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

extrait.leaf <- function(i,mask,image.fond.noir) {
  b <- bounding.rectangle(mask,i)
  leaf <- image.fond.noir[b$y,b$x,]
  mask.leaf <- mask[b$y,b$x]
  leaf[mask.leaf!=i] <- 0
  list(b=b,leaf=leaf)
}

## analyse d'une leaf
analyseLeaf <- function(x,lda1, lesion) {
  f <- x$leaf
  df6 <- data.frame(red=as.numeric(imageData(f)[,,1]), green=as.numeric(imageData(f)[,,2]), blue=as.numeric(imageData(f)[,,3]))
  df6$predict <- predict(lda1, df6)$class
  df6$tache <- as.numeric(df6$predict==lesion)
  df6$tache[df6$red+df6$green+df6$blue==0] <- 0
  mask <- channel(f, "gray")
  tache <- matrix(df6$tache, nrow=nrow(imageData(mask)))
  imageData(mask) <- tache

  ## dilatation
  brush <- makeBrush(lesionBorderSize, shape='disc')
  ##    mask <- dilateGreyScale(mask, brush)
  mask <- dilate(mask, brush)

  ## remplissage vides
  mask <- fillHull(mask)

  ## erosion
  mask <- erode(mask, brush)

  ## segmentation
  mask[mask<0] <- 0
  mask <- bwlabel(mask)
  features <- computeFeatures.shape(mask)

  ## suppression des petits objets
  w.petit <- which(features[,"s.area"]<lesionMinSize)
  mask <- rmObjects(mask,w.petit)

  features <- computeFeatures.shape(mask)
  list(features=features, mask=mask)
}

analyseFiles <- function(fileRdata=NA,pathResult,pathImages,onefileImage=NA,leafMinSize,leafBorderSize,lesionBorderSize,lesionMinSize,colorLesion) {
  ## chargement du résultat de l'analyse discriminante
  if (!is.na(fileRdata)) load(fileRdata)
  
  # return(list(fileRdata,pathResult,pathImages,onefileImage,leafMinSize,leafBorderSize,lesionBorderSize,lesionMinSize,colorLesion))
  # exit
  #ASSIGNE TO GLOBAL ENV
  fileRdata <<- fileRdata
  pathResult <<- pathResult
  pathImages <<- pathImages
  onefileImage <<- onefileImage
  leafMinSize <<- leafMinSize
  leafBorderSize <<- leafBorderSize
  lesionBorderSize <<- lesionBorderSize
  lesionMinSize <<- lesionMinSize
  colorLesion <<- colorLesion


  if (is.na(onefileImage)) onefileImage <- list.files(pathImages)
  nbfiles <- length(onefileImage)
  c <- 1
  for (image in onefileImage){
    incProgress(c/nbfiles, detail = paste("analysis leaf ", c, "/", nbfiles))
    analyseUniqueFile(pathResult,pathImages,image,leafMinSize,leafBorderSize,lesionBorderSize,lesionMinSize,colorLesion)
    c <- c + 1
  }
}

analyseUniqueFile <- function(pathResult,pathImages,onefileImage,leafMinSize,leafBorderSize,lesionBorderSize,lesionMinSize,colorLesion) {

  background <- names(lda1$prior)[1]
  limb <- names(lda1$prior)[2]
  lesion <- names(lda1$prior)[3]

  ## lecture de l'image source
  source.image <- paste(pathImages,'/',onefileImage,sep='')
  image <- readImage(source.image)
  ## imageData(image) <- imageData(image)[,,1:3]

  ## prédiction sur l'image floutée
  df5 <- data.frame(red=as.numeric(imageData(image)[,,1]), green=as.numeric(imageData(image)[,,2]), blue=as.numeric(imageData(image)[,,3]))
  df5$predict <- predict(lda1, df5)$class

  ## création des identificateurs des taches et des leafs
  df5$tache <- as.numeric(df5$predict==lesion)
  df5$leaf <- as.numeric(df5$predict!=background)

  ## masque des leafs
  mask <- channel(image, "gray")
  leaf <- matrix(df5$leaf, nrow=nrow(imageData(mask)))
  imageData(mask) <- leaf

  ## segmentation
  mask <- bwlabel(mask)
  features <- computeFeatures.shape(mask)

  ## suppression des objets plus petits que la surface minimum d'une leaf
  w <- which(features[,"s.area"]<leafMinSize)
  ## ligne suivante remplacée car il semble que les objets sont renumérotés dans certaines versions de R
  ## mask <- rmObjects(mask,w)
  mask[mask %in% w] <- 0

  ## suppression des vides
  mask <- fillHull(mask)

  ## suppression de la bordure par érosion
  brush <- makeBrush(leafBorderSize, shape='disc')
  mask <- erode(mask, brush)

  ## suppression du fond
  image.fond.noir <- image
  image.fond.noir[mask==0] <- 0

  ## séparation des leafs
  ## ligne suivante supprimée
  ## mask <- bwlabel(mask)
  features <- computeFeatures.shape(mask)
  ## ligne suivante remplacée car le numéro de l'objet ne correspond plus au numéro de ligne
  ## li <- lapply(1:nrow(features),extrait.leaf,mask,image.fond.noir)
  li <- lapply(as.numeric(row.names(features)),extrait.leaf,mask,image.fond.noir)

  ## analyse des leafs
  analyse.li <- lapply(li, analyseLeaf, lda1, lesion)

  filename <- strsplit(onefileImage,'\\.')[[1]][1]
  pdfname <- paste(filename,".pdf",sep='')
  txtname1 <- paste(filename,"_1.txt",sep='')
  txtname2 <- paste(filename,"_2.txt",sep='')
  if (!file.exists(pathResult)) dir.create(pathResult)
  pdffile <- paste(pathResult,'/',pdfname,sep='')
  txtfile1 <- paste(pathResult,'/',txtname1,sep='')
  txtfile2 <- paste(pathResult,'/',txtname2,sep='')

  pdf(pdffile)
  display(image, method="raster")

  ## sortie des résultats et coloration des lésions
  result <- NULL
  for (i in 1:length(li)) {
    result <- rbind(result,data.frame(fichier=filename,leaf=i,surfaceLeaf=features[i,"s.area"],surfaceLesion=if (is.null(analyse.li[[i]]$features)) 0 else analyse.li[[i]]$features[,"s.area"]))
    ## la ligne suivante a été remplacée par 3 lignes suite à une erreur apparue sur certaines versions de R
    ## image[li[[i]]$b$y,li[[i]]$b$x,][analyse.li[[i]]$mask>0] <- colorLesion
    tmpimage <- image[li[[i]]$b$y,li[[i]]$b$x,]
    tmpimage[analyse.li[[i]]$mask>0] <- colorLesion
    image[li[[i]]$b$y,li[[i]]$b$x,] <- tmpimage
  }
  row.names(result) <- NULL

  display(image, method="raster")
  dev.off()

  write.table(result,file=txtfile1,quote=FALSE,row.names=FALSE,sep='\t')

  ag.count <- aggregate(result$surfaceLesion,result[c("fichier","leaf", "surfaceLeaf")],length)
  names(ag.count)[4] <- "nb.lesions"
  ag.surface <- aggregate(result$surfaceLesion,result[c("fichier","leaf", "surfaceLeaf")],sum)
  names(ag.surface)[4] <- "surfaceLesions"
  ag <- merge(ag.count,ag.surface)
  ag$pourcent.lesions <- ag$surfaceLesions/ag$surfaceLeaf*100
  ag$nb.lesions[ag$surfaceLesions==0] <- 0

  write.table(ag[order(ag$leaf),],file=txtfile2,quote=FALSE,row.names=FALSE,sep='\t')
}
