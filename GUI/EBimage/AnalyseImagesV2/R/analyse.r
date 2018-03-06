## détection de lésions sur image couleur
## phase 2 : analyse d'image

source("fonctions_analyse.r")

## -------------------- Paramètres de l'analyse -----------------------------------
surface.feuille.mini <- 1000 ## surface minimum d'une feuille
bordure.feuille <- 3 ## épaisseur de bordure de feuille à supprimer
bordure.lesion <- 3 ## épaisseur de bordure de lésion à dilater / éroder
surface.lesion.mini <- 10 ## surface minimum d'une lésion
couleur.lesion <-  0 ## couleur des lésions dans l'image analysée (0=noir, 1=blanc)

## -------------------- Répertoires et fichiers -----------------------------------
path.image <- "../Images" ## Répertoire de stockage des fichiers images sources
file.image <- "IMG_5589_50.JPG" ## Fichier image source

path.sample <- "../Samples/5589" ## Répertoire de stockage des fichiers échantillons
background <- "fond" ## sous-répertoire contenant les échantillons de fond
limb <- "limbe" ## sous-répertoire contenant les échantillons de limbe
lesion <- "lesion" ## sous-répertoire contenant les échantillons de lesion

path.result <- "../Result" ## Répertoire de stockage des résultats d'analyse, créé si inexistant (peut être le même que path.image)

## --------------------- Analyse --------------------------------------------------
source("analyse_image.r")

## ------------- Fin d'analyse ----------------------------------------------------

## ------------- Les lignes suivantes sont des exemples d'autres analyses ---------
path.sample <- "../Samples/5583" ## Répertoire de stockage des fichiers échantillons
file.image <- "IMG_5583_50.JPG" ## Fichier image source

path.sample <- "../Samples/5605" ## Répertoire de stockage des fichiers échantillons
file.image <- "IMG_5605_50.JPG" ## Fichier image source

path.sample <- "../Samples/5605" ## Répertoire de stockage des fichiers échantillons
file.image <- "IMG_5609_50.JPG" ## Fichier image source

path.sample <- "../Samples/5605" ## Répertoire de stockage des fichiers échantillons
file.image <- "IMG_5598_50.JPG" ## Fichier image source

path.sample <- "../Samples/5605" ## Répertoire de stockage des fichiers échantillons
file.image <- "IMG_5581_50.JPG" ## Fichier image source

path.sample <- "../Samples/5605" ## Répertoire de stockage des fichiers échantillons
file.image <- "IMG_5583_50.JPG" ## Fichier image source

path.sample <- "../Samples/5605" ## Répertoire de stockage des fichiers échantillons
file.image <- "IMG_5593_50.JPG" ## Fichier image source

path.sample <- "../Samples/5605" ## Répertoire de stockage des fichiers échantillons
file.image <- "IMG_5591_50.JPG" ## Fichier image source

path.sample <- "../Samples/5615" ## Répertoire de stockage des fichiers échantillons
file.image <- "IMG_5615_50.JPG" ## Fichier image source

## ----------------- Fin de fichier -----------------------------------------------
