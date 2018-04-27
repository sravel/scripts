## détection de lésions sur image couleur
## phase 2 : analyse d'image

source("fonctions_analyse.r")

## -------------------- Paramètres de l'analyse -----------------------------------
surface.feuille.mini <- 1000 ## surface minimum d'une feuille
bordure.feuille <- 3 ## épaisseur de bordure de feuille à supprimer
bordure.lesion <- 3 ## épaisseur de bordure de lésion à dilater / éroder
surface.lesion.mini <- 10 ## surface minimum d'une lésion
couleur.lesion <-  0 ## couleur des lésions dans l'image analysée (0=noir, 1=blanc)

## -------------------- Répertoires et fichiers Exemple1---------------------------
fileRdata <- "../Exemple1/Samples/Samples.RData" ## Répertoire de stockage des fichiers échantillons
path.result <- "../Exemple1/Result"  ## Répertoire de stockage des résultats d'analyse, créé si inexistant (peut être le même que path.image)
path.image <- "../Exemple1/Images"   ## Répertoire de stockage des fichiers images sources
file.image <- "IMG_5593_50.jpg"      ## Fichier image source

## -------------------- Répertoires et fichiers Exemple2 --------------------------
fileRdata <- "/media/sebastien/Bayer/AnalyseImagesV4/Samples/Samples.RData"   ## Répertoire de stockage des fichiers échantillons
path.result <- "/media/sebastien/Bayer/AnalyseImagesV4/Exemple2/Result"    ## Répertoire de stockage des résultats d'analyse, créé si inexistant (peut être le même que path.image)
path.image <- "/media/sebastien/Bayer/AnalyseImagesV4/Images"     ## Répertoire de stockage des fichiers images sources


## -------- Exemple analyse d'un répertoire complet -------------------------------
analyse.image(fileRdata=fileRdata,
              path.result=path.result,
              path.image=path.image,
              file.image=NA, ## analyse du répertoire complet
              surface.feuille.mini=surface.feuille.mini,
              bordure.feuille=bordure.feuille,
              bordure.lesion=bordure.lesion,
              surface.lesion.mini=surface.lesion.mini,
              couleur.lesion=couleur.lesion)
## ------------- Fin d'analyse ----------------------------------------------------

## ----------------- Fin de fichier -----------------------------------------------
