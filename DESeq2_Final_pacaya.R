#source: Analyse statistique des données RNA Seq, DESeq2.  J. Aubert, A. de la Foye et C. Hennequet-Antier INRA 2013 
#Updated by: Abdoulaye DIALLO, Christine TRANCHANT, Sebastien RAVEL 

###################################################
###code1:InstalLibrairies
###################################################
#source("http://www.bioconductor.org/biocLite.R")
#biocLite(pkgs=c("DESeq","DESeq2","Biobase"))

# Chargement des librairies
library(DESeq2)
library(Biobase)


# Déplacement dans le répertoire de travail
setwd("/home/abdou/Documents/Stage_pacaya_2016/R/analyse_diff_DMN6S_WMN6C_g33_par_g10/")

# Prefixe qui sera ajouté à chaque image dans le repertoire 
prefix=""

###################################################
### code2:Chargement de la table de comptage
###################################################

# creation dun object countTable  contenant les counts obtenus à partir des fichiers HTSEQCOUNT
countTable <- read.table("../../g02L5-6_htseqCount.csv", header=TRUE, row.names=1, sep = "\t")
head(countTable)
#colnames(countTable)

## On selectionne les colonnes qui nous interessent: les individus correspondant aux codes DMN6S et WMN6C
#colonnes_DMN6S <- c("g09L005","g09L006","g21L005","g21L006","g33L005","g33L006")
#colonnes_WMN6C <- c("g38L005","g38L006","g40L005","g40L006","c03L005","c03L006")

# Pour eviter les billets qui peuvent être liés aux reliquats de nature technique, on fait un choix alterné
# des indivdus sauvages et cultivés par lane.
colonnes_L5 <- c("g09L005","g38L005","g21L005")
colonnes_L6 <- c("g40L006","g10L006","c03L006")

head(countTable[,c(colonnes_L5,colonnes_L6)])
countTable <- countTable[,c(colonnes_L5,colonnes_L6)]

###################################################
### Code3: Construction objet DESeqDataSet:dds
###################################################
##création de la table colData qui décrit le plan d’expérience:
##"Liste des conditions comparées"
## les facteurs sont des structures de données semblables auc vectors. type de données hétérogènes.
#condition<-factor(c(rep("wild",8),rep("domesticated",10),rep("domesticated",12),rep("wild",2),rep("domesticated",2),rep("wild",4)))
condition<-factor(rep(c("domesticated","wild"),3))

##Type: single-read ou paired-end
type<-factor(rep("paired-end",6))        # rep("A",x): repetition de x fois A

## Création de la data.frame (tableaux de données): matrice de données de type differents
colData<-data.frame(condition,type,row.names=colnames(countTable)) 
head(colData)

##création de l’objet dds
dds<-DESeqDataSetFromMatrix(countTable,colData,design=~condition)

#objetdds
class(dds)
colData(dds)
design(dds)

#manipulation des données de comptage
# vérifier que le nombre de colonnes et de lignes correspondent bien aux données attendues.
dim(counts(dds))
head(counts(dds))
summary(counts(dds))

###################################################
### code4: Exploration des données
###################################################

##combien il y a−t−il de comptages nuls par échantillon?
apply(counts(dds),2,FUN=function(x)sum(x==0))

## combien de comptages nuls en % ?
apply(counts(dds), 2, FUN = function(x) sum(x ==0))/nrow(counts(dds))

##quelles sont les profondeurs de séquençage de chaque échantillon?
colSums(counts(dds))



#png(paste(prefix,"profondeurs de séquençage de chaque échantillon.png"), width = 1500, height = 1500, res=200)
#barplot(colSums(counts(dds)))
#dev.off()

##distribution des gènes suivant leur log−fréquence pour chaque échantillon
png(paste(prefix,"4a_frequenceDeComptage_avantFiltre.png"), width = 1500, height = 1500, res=200)
par(mfrow=c(2,3))  # constitue un documents de 6 images (2 lignes, 3 colonnes)
for(i in 1:6){
  hist(log(counts(dds)[,i]+1),
       breaks=seq(0,20,1),col="grey",cex.lab=0.8,
       main=colnames(counts(dds))[i],
       xlab="Valeur du comptage (nb reads/gènes) log(count+1)",
       ylab="Fréquence du comptage")
}
dev.off()

###################################################
### code5:Filtrage
###################################################
## filtrage arbitraire: on ne garde que les gènes avec un comptage moyen supérieur à 5 (seuil minumum de couverture).
count_mean5<-countTable[rowMeans(countTable)>=5,]
#count_sum5<-countTable[rowSums(countTable)>=5,]
#dds_sum<-DESeqDataSetFromMatrix(count_sum5,colData,design=~condition)
#dim(counts(dds_sum))
dds<-DESeqDataSetFromMatrix(count_mean5,colData,design=~condition)
dim(counts(dds))

png(paste(prefix,"5a_frequenceDeComptage_apresFiltre_5_reads.png"), width = 1500, height = 1500, res=200)
par(mfrow=c(2,3))
for(i in 1:6){
  hist(log(counts(dds)[,i]+1),
       breaks=seq(0,20,1),col="grey",cex.lab=0.8,
       main=colnames(counts(dds))[i],
       xlab="Valeur du comptage (nb reads/gènes) log(count+1)",
       ylab="Fréquence du comptage")
}
dev.off()

head(counts(dds))

#enregistre dans un tableau les comptage filtre
write.table(count_mean5, paste(prefix,"5b_tableComptage_filtre_5_reads.csv"), sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

###ACP
## Permet de regarder si l'expression des genes est très differente ou si les replicats sont suffisament semblables. Dans le cas contraire eliminer les mauvais replicats.
par(mfrow=c(1,2))
acp <- prcomp(t(count_mean5))            # Analyse de la composante principale de la transposée(echange de lignes et de colonnes) sur les données filtrées.
png(paste(prefix,"5c_ACP_filtre_5_reads_biplot1.png"), width = 1500, height = 1500, res=200)
biplot(acp, var.axes=FALSE)              # Tracée d'un biplot (plot multidimensionel)
dev.off()

png(paste(prefix,"5e_ACP_filtre_5_reads_biplot2.png"), width = 1500, height = 1500, res=200)
biplot(acp,choices = 1:2, scale = 1, pc.biplot = FALSE)
dev.off()

#####################################################
### code6: Normalisation et Analyse différentielle
#####################################################
## la fonction DESeq() réalise la normalisation et l’analyse diff en une seule étape
## la normalisation des données de comptage se fait par la méthode "effective library size".
dds <- DESeq(dds)               # Analyse de l'expression differentielle basée sur la distribution négative binomiale
colData(dds)
#DESeq2 calcule un size factor Sj par échantillon; afin de normaliser les comptages: 
# le nb de reads pour le gene i dans l'echantillon j noté X'ij = Xij/Sj
sizeFactors(dds)
# Hypothèses:
# La majorité des genes n'est pas differentiellement exprimées.
#Autant de gènes sur-exprimés que de sous-exprimés.

### code6.1: Visualisation de la normalisation ###

## quelles sont les profondeurs de séquençage de chaque échantillon?
colSums(counts(dds,normalized=FALSE))

## profondeur de séquençage après normalisation
colSums(counts(dds,normalized=TRUE))



## Visualisation graphique Avant/Après

png(paste(prefix,"6a_profondeurSeqFiltre5_vs_Normalised.png"), width = 3000, height = 1500, res=200)
par(mfrow=c(1,2))      # mettre par(mfrow=c(1,1)) pour 2 graphiques
barplot(colSums(counts(dds)),
        main="profondeur de séquençage avant normalisation")
barplot(colSums(counts(dds,normalized=TRUE)),
        main="profondeur de séquençage après normalisation")
dev.off()

#données de comptage normalisées
head(counts(dds,normalized=TRUE),n=6)
summary(counts(dds,normalized=TRUE))

##distribution des gènes suivant leur log−fréquence pour chaque échantillon après normalisation
png(paste(prefix,"6b_frequenceDeComptageFiltre5_post_normalized.png"), width = 3000, height = 1500, res=200)
par(mfrow=c(2,3))
for(i in 1:6){
  hist(log(counts(dds,normalized=TRUE)[,i]+1),
       breaks=seq(0,20,1),col="grey",cex.lab=0.8,
       main=colnames(counts(dds))[i],
       xlab="Valeur du comptage normalisé (nb reads/gènes) log(count+1)",
       ylab="Fréquence du comptage")
}
dev.off()

# manipulation des données de comptage
dim(counts(dds,normalized=TRUE))
head(counts(dds,normalized=TRUE))
print(summary(counts(dds,normalized=TRUE)))


### code6.2: Dispersion des gènes###

#La fonction DESeq() a effectué une estimation de la dispersion des gènes
#Valeurs de dispersion:
disp<-as.data.frame(dispersions(dds))
head(disp,n=3)
print(summary(disp))

#Boxplot des valeurs de dispersion
## The dispersion can be understood as the square of the coefficient of biological variation. 
#So, if a gene’s expression typically differs from replicate to replicate sample by 20%, 
#this gene’s dispersion is .202 = .04. Note that the variance seen between counts is the sum of 
#two components: the sample-to-sample variation just mentioned, and the uncertainty in measuring a 
#concentration by counting reads
png(paste(prefix,"6c_dispersionGene.png"), width = 1500, height = 1500, res=200)
par(mfrow=c(1,1))
boxplot(sqrt(disp),horizontal=T,las=1,cex=0.5)
title("racine carrée de la dispersion calculée par DESeq2")
dev.off()

### code6.3: relation entre la dispersion et la moyenne des comptages ###
png(paste(prefix,"6d_dispertionMeanCount.png"), width = 3000, height = 1500, res=200)
par(mfrow=c(1,1))
DESeq2::plotDispEsts(dds)       # estimation de la dispersion
title("relation entre la dispersion et la moyenne des comptages")
dev.off()

### code6.4: details sur les calculs des valeurs de dispersion, tests statistiques et pvalues ###
## p-value p: probabilité d'obtenir la même valeur (ou une valeur encore plus extrême) 
## du test si l'hypothèse nulle était vraie.  En général, les hypothèses s'opposant à l'hypothèse nulle ont le fardeau de la preuve
## Dans notre cas, l'hypothèse nulle serait qu'il n'y ait pas de difference de genes entre domestiqués et sauvages.
## p <= 0.01 : très forte présomption contre l'hypothèse nulle
# 0.01 < p ≤ 0.05 : forte présomption contre l'hypothèse nulle
# 0.05 < p ≤ 0.1 : faible présomption contre l'hypothèse nulle
# p > 0.1 : pas de présomption contre l'hypothèse nulle
details<-mcols(dds,use.names=T)
head(details,n=3)
mcols(mcols(dds,use.names=T))

###code6.5: Graphique MA−plot: relation entre le comptage moyen d’un gène et son log2−ratio entre les 2 conditions ###
png(paste(prefix,"6c_MAplot.png"), width = 2000, height = 2000, res=200)
par(mfrow=c(1,1))
DESeq2::plotMA(dds, main="DESeq2", ylim=c(-2,2))
dev.off()

# les gènes différentiellement exprimés au seuil de 10% après ajustement
# des tests multiples par procédure de Benjamini−Hochberg(BH) sont indiqués
# en rouge

#####################################################
### code7: Résultats de l’analyse différentielle
#####################################################
#Pour chaque gène i, DESeq2 et edgeR donnent : 
#  une estimation de log2(FCi), donc de FCi
#  la précision de cette estimation (écart-type)
# donc la p-valeur associée au gène i
#L’ensemble des N  p-valeurs obtenues est ajusté afin de conclure

# res contient les résultats de l’analyse différentielle pour la dernière variable
res<-results(dds)
head(res)
resultsNames(dds)
mcols(res,use.names=T)
summary(res)

resOrdered <- res[order(res$padj),]             # on s'interesse à la p-value
head(resOrdered)

# distribution des p−valeurs brutes
png(paste(prefix,"7a_p-valeursBrut.png"), width = 2000, height = 1500, res=200)
par(mfrow=c(1,1))
hist(res$pvalue,breaks=100,border="slateblue", col ="orange",
     main="Histogramme des p−values brutes")
dev.off()

## Avec les P-value ajustées
# Combien de gènes sont déclarés DE à 5%?
table(res$padj<0.05,useNA="always")

# Selection des gènes déclarés DE à 5% dans une table ordonnée par p−value ajustées croissante
resSig5<-na.omit(res)
resSig5<-DESeq2::results(dds, alpha=0.05)
resSig5<-resSig5[order(resSig5$padj),]
nrow(resSig5)
head(resSig5)

# Combien de gènes sont déclarés DE à 1%?
table(res$padj<0.01,useNA="always")

# Selection des gènes déclarés DE à 1% dans une table ordonnée par p−value ajustées croissante
resSig1<-na.omit(res)
resSig1<-DESeq2::results(dds, alpha=0.01)
resSig1<-resSig1[order(resSig1$padj),]
nrow(resSig1)
head(resSig1)

# selection des genes significatifs  au seuil alpha=0.01 et .05
deseq2_resSig <- na.omit(res)
deseq2_resSig_5p <- res[which(res$padj < 0.05),]
deseq2_resSig <- res[which(res$padj < 0.01),]
head(deseq2_resSig_5p)
head(deseq2_resSig)


##enregistrement dans un tableau
write.csv(as.data.frame(resOrdered), file=paste(prefix,"7b_condition_domesticated_results"), sep="\t")
write.table (deseq2_resSig,
             file=paste(prefix,"7c_domesticated_vs_wild_0.01"),
             quote=FALSE, row.names=TRUE, sep="\t")

write.table (deseq2_resSig_5p,
             file=paste(prefix,"7d_domesticated_vs_wild_0.05"),

             quote=FALSE, row.names=TRUE,  sep="\t")

table(res$padj < 0.05,useNA="always")
table(res$padj < 0.01,useNA="always")

# sep Up et Down
# fold Change mesure combien de fois une quantité a changé entre 2 conditions 
#(fold = nombre de fois, ex: valeur initiale à A et final à B ==> fold change = (B-A)/A.)
deseq2_resSig_5pUp <- deseq2_resSig_5p[which(deseq2_resSig_5p$log2FoldChange > 1),]
nrow(deseq2_resSig_5pUp)
write.table (deseq2_resSig_5pUp,
             file=paste(prefix,"7e_UP_0.05"),quote=FALSE, row.names=TRUE, sep="\t")

deseq2_resSig_5pDown <- deseq2_resSig_5p[which(deseq2_resSig_5p$log2FoldChange < -1),]
nrow(deseq2_resSig_5pDown)
write.table (deseq2_resSig_5pDown,
             file=paste(prefix,"7f_Down_0.05"),quote=FALSE, row.names=TRUE, sep="\t")


#######################################################################
### code8: Comparaison des p−values brutes, ajustées Bonferroni et BH
#######################################################################
#calcul des p−values ajustées selon la méthode de bonferroni
res$padjBonf<-p.adjust(res$pval,meth="bonferroni")
head(res)

#histogramme des distributions des p−values
png(paste(prefix,"8a_distribution-p-values.png"), width = 2000, height = 1500, res=200)
par(mfrow=c(1,1))
hist(res$padjBonf,border="blue",xlab="pvalue",
     main="distributions␣des␣p−values")
hist(res$padj,add=TRUE,border="red")
hist(res$pval,add=TRUE)
legend(x="topleft",legend=c("pVal","padjBonf","padjBH"),
       fill=c("black","blue","red"))
dev.off()


#cumul des p−values
png(paste(prefix,"8b_distributionCumul-p-values.png"), width = 2000, height = 1500, res=200)
resord<-res[order(res$pval),]
head(resord)
plot(resord$pval,ylab="Pvalues",cex=0.5,xlab="rang␣du␣gène")
points(resord$padjBonf,col="blue",cex=0.5)
points(resord$padj,col="red",cex=0.5)
abline(h=0.05,lty=2)
legend(x="bottomright",legend=c("pVal","padjBonf","padjBH"),
       col=c("black","blue","red"),pch=1)
dev.off()


###########################################################
###code9:ClassificationAscendanteHiérarchique(CAH)
###########################################################
#CAH
cdslog <- log (counts(dds)+1)
dist.cor <- 1-cor(cdslog , use="pairwise.complete.obs", method = "spearman")
hc.cor <- hclust (as.dist (dist.cor), method="ward.D")

#graph Phylogénie
png(paste(prefix,"9a_CAH_phylo.png"), width = 1500, height = 1500, res=200)
plot (hc.cor, main=paste("CAH DES echantillons à partir des comptages \n de ",
                         nrow(cdslog), " genes ", sep=""),
      sub="DIstance= 1-correlation")
dev.off()

colData(dds)



###########################################################
###code9bis ?????
###########################################################

#visualisation et controles qualité
rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

rlogMat <- assay(rld)
vstMat <- assay(vsd)

png(paste(prefix,"9b_DESeq2_VST_and_log2.png"), width = 1500, height = 1500, res=200)
par(mfrow=c(1,1))
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord] < 150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c("blue", "black")
matplot(px[ord], cbind(assay(vsd)[, 1], log2(px))[ord, ], type="l", lty=1, col=vstcol, xlab="n", ylab="f(n)")
legend("bottomright", legend = c(expression("variance stabilizing transformation"), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()
#dev.copy(png,"DESeq2_VST_and_log2.png")


###########################################################
### code10: HEATMAP
###########################################################

library (gplots)
library(RColorBrewer)

#Stabilization variance is usefull for visualization methods like clustering,PCA
rld<-rlogTransformation(dds,blind=T)
dists<-dist(t(assay(rld)))
mat<-as.matrix(dists)
hmcol=colorRampPalette(brewer.pal(9,"GnBu"))(100)


png(paste(prefix,"10a_heatmap1.png"), width = 4000, height = 4000, res=300)
heatmap.2(mat,trace="none",col=rev(hmcol),margin=c(10,6),
          main="heatmap with stabilization variance regularized log Transformation")
dev.off()

# Selection du nombre de gène à afficher (10)
#select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

select <- order(res$pval)[1:20] #for 1st 20 genes
select <- order(res$pval)[1:621] #for 1st 20 genes
select <- order(res$pval)[1:1109] #for 1st 20 genes

# HEATMAP
par(mfrow=c(1,1))


png(paste(prefix,"10b_heatmap2.png"), width = 8000, height = 10000, res=500)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(5,10))
dev.off()


png(paste(prefix,"10c_heatmap3.png"), width = 4000, height = 5000, res=300)

heatmap.2(counts(dds,normalized=TRUE)[select,], col=redgreen(75), scale="row",  key=TRUE,dendrogram="both",
          y=FALSE,density.info="none", cexRow=0.6,cexCol=1,margins=c(5,10),  trace="none",srtCol=45)
dev.off()


png(paste(prefix,"10d_heatmap4.png"), width = 4000, height = 5000, res=300)
heatmap.2(rlogMat[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(5, 10))
dev.off()


png(paste(prefix,"10e_heatmap5.png"), width = 4000, height = 5000, res=300)
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(5, 10))
dev.off()





##### TEST MGG_11605T0
