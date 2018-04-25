## détection de lésions sur image couleur
## phase 1 : apprentissage à partir d'images échantillonnées
## path.sample est le répertoire contenant les sous-répertoires contenant les fichiers d'échanntillons
## les trois derniers arguments de la fonction apprentissage sont (dans cer ordre)
## chaque sous-réopertoire contient un nombre indéterminé de fichiers d'une même catégorie (fond ou limbe ou lésion)
## trois sous-répertoires contenant trois catégories de pixels (fond, limbe, lésion) sont requis,
## répertoire path.sample peut contenir d'autres sous-répertoires inutilisés

source("fonctions_apprentissage_V4.r")

## Choisir Exemple1 ou Exemple2
path.sample <- "../Exemple1/Samples"
path.sample <- "../Exemple2/Samples"

apprentissage(path.sample,"fond","limbe","lesion")

## Fin de fichier
