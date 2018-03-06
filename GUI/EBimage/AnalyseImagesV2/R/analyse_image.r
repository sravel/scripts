## chargement du résultat de l'analyse discriminante
filename <- paste0(tail(strsplit(path.sample,'/')[[1]],1),".RData")
file.lda <- paste(path.sample,'/',filename,sep='')
load(file.lda)
lda1

## lecture de l'image source
source.image <- paste(path.image,'/',file.image,sep='')
image <- readImage(source.image)
## imageData(image) <- imageData(image)[,,1:3]

## prédiction sur l'image floutée
df5 <- data.frame(red=as.numeric(imageData(image)[,,1]), green=as.numeric(imageData(image)[,,2]), blue=as.numeric(imageData(image)[,,3]))
df5$predict <- predict(lda1, df5)$class

## création des identificateurs des taches et des feuilles
df5$tache <- as.numeric(df5$predict==lesion)
df5$feuille <- as.numeric(df5$predict!=background)

## masque des feuilles
mask <- channel(image, "gray")
feuille <- matrix(df5$feuille, nrow=nrow(imageData(mask)))
imageData(mask) <- feuille

## segmentation
mask <- bwlabel(mask)
features <- computeFeatures.shape(mask)

## suppression des objets plus petits que la surface minimum d'une feuille
w <- which(features[,"s.area"]<surface.feuille.mini)
mask <- rmObjects(mask,w)

## suppression des vides
mask <- fillHull(mask)

## suppression de la bordure par érosion
brush <- makeBrush(bordure.feuille, shape='disc')
mask <- erode(mask, brush)

## suppression du fond
image.fond.noir <- image
image.fond.noir[mask==0] <- 0

## séparation des feuilles
features <- computeFeatures.shape(mask)
li <- lapply(1:nrow(features),extrait.feuille)

## analyse des feuilles
analyse.li <- lapply(li, analyse.feuille)

filename <- strsplit(file.image,'\\.')[[1]][1]
pdfname <- paste(filename,".pdf",sep='')
txtname1 <- paste(filename,"_1.txt",sep='')
txtname2 <- paste(filename,"_2.txt",sep='')
if (!file.exists(path.result)) dir.create(path.result)
pdffile <- paste(path.result,'/',pdfname,sep='')
txtfile1 <- paste(path.result,'/',txtname1,sep='')
txtfile2 <- paste(path.result,'/',txtname2,sep='')

pdf(pdffile)
display(image, method="raster")

## sortie des résultats et coloration des lésions
result <- NULL
for (i in 1:length(li)) {
    result <- rbind(result,data.frame(fichier=filename,feuille=i,surface.feuille=features[i,"s.area"],surface.lesion=if (is.null(analyse.li[[i]]$features)) 0 else analyse.li[[i]]$features[,"s.area"]))
    imageData(image[li[[i]]$b$y,li[[i]]$b$x,])[imageData(analyse.li[[i]]$mask)>0] <- couleur.lesion
}
row.names(result) <- NULL

display(image, method="raster")
dev.off()

write.table(result,file=txtfile1,quote=FALSE,row.names=FALSE,sep='\t')

ag.count <- aggregate(result$surface.lesion,result[c("fichier","feuille", "surface.feuille")],length)
names(ag.count)[4] <- "nb.lesions"
ag.surface <- aggregate(result$surface.lesion,result[c("fichier","feuille", "surface.feuille")],sum)
names(ag.surface)[4] <- "surface.lesions"
ag <- merge(ag.count,ag.surface)
ag$pourcent.lesions <- ag$surface.lesions/ag$surface.feuille*100
ag$nb.lesions[ag$surface.lesions==0] <- 0

write.table(ag[order(ag$feuille),],file=txtfile2,quote=FALSE,row.names=FALSE,sep='\t')
