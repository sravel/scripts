options(width=160)

## packages nécessaires
library(EBImage)
library(MASS)

## retourne les indices extrêmes de la valeur objet dans le vecteur x
## (appelé par bounding.rectangle)
range.na <- function(x, object) {
    w <- which(x==object)
    if (length(w)==0) return(c(NA,NA))
    return(range(w))
}

bounding.rectangle <- function(mask, object) {
    m <- imageData(mask)
    range.x <- range(apply(m,1,range.na,object),na.rm=TRUE)
    range.y <- range(apply(m,2,range.na,object),na.rm=TRUE)
    list(x=range.x[1]:range.x[2],y=range.y[1]:range.y[2])
}

extrait.feuille <- function(i) {
    b <- bounding.rectangle(mask,i)
    feuille <- image.fond.noir[b$y,b$x,]
    mask.feuille <- mask[b$y,b$x]
    feuille[mask.feuille!=i] <- 0
    list(b=b,feuille=feuille)
}

## analyse d'une feuille
analyse.feuille <- function(x) {
    f <- x$feuille
    df6 <- data.frame(red=as.numeric(imageData(f)[,,1]), green=as.numeric(imageData(f)[,,2]), blue=as.numeric(imageData(f)[,,3]))
    df6$predict <- predict(lda1, df6)$class
    df6$tache <- as.numeric(df6$predict==lesion)
    df6$tache[df6$red+df6$green+df6$blue==0] <- 0
    mask <- channel(f, "gray")
    tache <- matrix(df6$tache, nrow=nrow(imageData(mask)))
    imageData(mask) <- tache

    ## dilatation
    brush <- makeBrush(bordure.lesion, shape='disc')
    mask <- dilateGreyScale(mask, brush)

    ## remplissage vides
    mask <- fillHull(mask)

    ## erosion
    mask <- erode(mask, brush)

    ## segmentation
    mask[mask<0] <- 0
    mask <- bwlabel(mask)
    features <- computeFeatures.shape(mask)

    ## suppression des petits objets
    w.petit <- which(features[,"s.area"]<surface.lesion.mini)
    mask <- rmObjects(mask,w.petit)

    features <- computeFeatures.shape(mask)
    list(features=features, mask=mask)
}

