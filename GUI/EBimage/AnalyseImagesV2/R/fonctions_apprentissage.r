options(width=160)

## packages nécessaires
library(EBImage)
library(lattice)
library(MASS)

## fonction de lecture des images d'un groupe ; retourne le data.frame des pixels
load.group <- function(g,path.sample) {
    path.group <- paste(path.sample,g,sep='/')
    files.group <- list.files(path.group,full.name=TRUE)
    sample <- lapply(files.group,readImage)
    ## constitution du data frame des pixels échantillonnés
    li <- lapply(sample, function(im) {
        data.frame(group=g,red=as.numeric(imageData(im)[,,1]), green=as.numeric(imageData(im)[,,2]), blue=as.numeric(imageData(im)[,,3]))
    })
    do.call(rbind, li)
}

apprentissage <- function(path.sample) {

    ## Recherche des sous-répertoires de path.sample
    group <- list.dirs(path.sample,full.names=FALSE)[-1] ## -1 pour supprimer le premier nom (toujouts vide)

    ## constitution du data.frame des pixels des échantillons
    li <- lapply(group,load.group,path.sample)
    df2 <- do.call(rbind, li)

    ## analyse discriminante
    lda1 <- lda(df2[2:4], df2$group)

    ## nom commun aux 3 fichiers de sortie, identique au nom du réprtoire
    filename <- tail(strsplit(path.sample,'/')[[1]],1)

    ## écriture du fichier texte des résultats
    file.txt <- paste(path.sample,paste0(filename,".txt"),sep='/') ## fichier de sortie texte
    sink(file.txt)
    print(table(df2$group))
    print(lda1$scaling)
    df2$predict <- predict(lda1, df2[2:4])$class
    print(table(df2$group, df2$predict))
    sink()

    ## graphe des groupes dans le plan discriminant
    file.pdf <- paste(path.sample,paste0(filename,".pdf"),sep='/') ## fichier de sortie pdf
    df4 <- cbind(df2, as.data.frame(as.matrix(df2[2:4])%*%lda1$scaling))
    pdf(file.pdf)
    print(xyplot(LD2~LD1, group=group, cex=0.8, alpha=1, pch=1, asp=1, auto.key=TRUE, data=df4))
    dev.off()

    ## sauvegarde de l'analyse
    file.lda <- paste(path.sample,paste0(filename,".RData"),sep='/')
    save(lda1,file=file.lda)
}

## Fin de fichier
