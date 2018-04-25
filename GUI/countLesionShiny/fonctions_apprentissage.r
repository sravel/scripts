options(width=160)

## packages nécessaires
library(EBImage)
library(lattice)
library(MASS)

## fonction de lecture des images d'un groupe ; retourne le data.frame des pixels
load_group <- function(g,pathCalibration) {
  path_group <- paste(pathCalibration,g,sep='/')
  files_group <- list.files(path_group,full.name=TRUE)
  sample <- lapply(files_group,readImage)
  ## constitution du data frame des pixels échantillonnés
  li <- lapply(sample, function(im) {
    data.frame(group=g,red=as.numeric(imageData(im)[,,1]), green=as.numeric(imageData(im)[,,2]), blue=as.numeric(imageData(im)[,,3]))
  })
  do.call(rbind, li)
}

apprentissage <- function(pathCalibration,...) {
  ## Les arguments passés dans "..." doivent être (dans cet ordre) le nom (relatif) des sous-répertoires fond, limbe, lésions
  ## Recherche des sous-répertoires de pathCalibration
  dirs <- list.dirs(pathCalibration,full.names=FALSE)[-1] ## -1 pour supprimer le premier nom (toujouts vide)
  
  ## vérification de l'existence des sous-répertoires passés en argument
  group <- list(...)
  if (any(is.na(match(unlist(group),dirs)))) stop("Répertoire(s) inexistant(s).")
  
  ## constitution du data.frame des pixels des échantillons
  li <- lapply(group,load_group,pathCalibration)
  df2 <- do.call(rbind, li)
  
  ## analyse discriminante
  lda1 <- lda(df2[2:4], df2$group)
  
  ## nom commun aux 3 fichiers de sortie, identique au nom du réprtoire
  basename <- tail(strsplit(pathCalibration,'/')[[1]],1)
  
  # incProgress(2/3, detail = "make files")

  ## écriture du fichier texte des résultats
  file.txt <- paste(pathCalibration,paste0(basename,".txt"),sep='/') ## fichier de sortie texte
  sink(file.txt)
  print(table(df2$group))
  print(lda1$scaling)
  df2$predict <- predict(lda1, df2[2:4])$class
  sink()
  
  outCalibrationCSV <<- paste(pathCalibration,paste0(basename,"_info.csv"),sep='/') ## fichier de sortie csv
  outCalibrationTable <<- as.data.frame.matrix(table(df2$group, df2$predict))
  write.csv2(outCalibrationTable, file = outCalibrationCSV)
  
  ## graphe des groupes dans le plan discriminant
  plotFileCalibration <<- paste(pathCalibration,paste0(basename,".png"),sep='/') ## fichier de sortie png
  df4 <- cbind(df2, as.data.frame(as.matrix(df2[2:4])%*%lda1$scaling))
  
  png(plotFileCalibration)
  print(xyplot(LD2~LD1, group=group, cex=0.8, alpha=1, pch=1, asp=1, auto.key=TRUE, data=df4))
  dev.off()
  
  ## sauvegarde de l'analyse
  fileRData <<- paste(pathCalibration,paste0(basename,".RData"),sep='/')
  save(lda1,file=fileRData)
  return(list(code = 1, mess = fileRData))
}

## Fin de fichier
