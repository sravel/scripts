## détection de lésions sur image couleur
## phase 1 : apprentissage à partir d'images échantillonnées
## path.sample est le répertoire contenant les sous-répertoires contenant les fichiers d'échanntillons
## chaque sous-réopertoire contient u nombre indéterminé de fichiers d'une même catégorie (fond ou limbe ou lésion)
## tous les sous-répertoires sont analysés, seuls les sous-répertoires utiles doivent figurer dans path.sample
## les noms des sous-répertoires sont libres et doivent être spécifiés lors de la phase 2 (analyse)

source("fonctions_apprentissage.r")

## Répertoire des sous-répertoires d'échantillons
path.sample <- "../Samples/5589"
apprentissage(path.sample)
## Fin d'apprentissage

## Autres exemples
path.sample <- "../Samples/5589"
path.sample <- "../Samples/5605"
path.sample <- "../Samples/5615"

## Fin de fichier
