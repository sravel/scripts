#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package run_DAPC.py
# @author Sebastien Ravel

"""
	The run_DAPC script
	===================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 08/07/2016
	:version: 0.1

	Script description
	------------------

	This Programme run DAPC for k i to m

	Example
	-------

	>>> run_DAPC.py -m asie_480mlg.txt -o asie_480mlg_DAPC.txt -pi 2 -pm 10 -pca 20

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display run_DAPC.py version number and exit

	Input mandatory infos for running:
		- \-m <filename>, --mat <filename>
						matrice file path
		- \-o <filename>, --order <filename>
						file with re-order name of matrice

	Input infos for running with default values:
		- \-pca <int>, --pcanum <int>
						Number value of PCA retains (default = NULL)
		- \-pi <int>, --popi <int>
						Number of pop Min (default = 2)
		- \-pm <int>, --popm <int>
						Number of pop Max (default = 10)

"""

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import replace_all, relativeToAbsolutePath, loadInListCol, loadInDictLine, existant_file

## Python modules
import argparse
from time import localtime, strftime
from shutil import rmtree
##################################################
## Variables Globales
version="0.3"
VERSION_DATE='21/04/2016'

##################################################
## PARAM DAPC R
DAPCfix = """

f1 <- function(vec) {
  m <- mean(vec, na.rm = TRUE)
  vec[is.na(vec)] <- m
  return(vec)
}


library(ade4)
library(adegenet)
library(ggplot2)

mvdata <- read.table("**PATHTOFILE**", header=TRUE,sep="\\t",na.string="-9")
attach(mvdata)

vect=c(rep(1,**INDIV**))
name= mvdata[c(1:**INDIV**),1]

W <- df2genind(X=mvdata[,c(2:**MAKERS**)],pop=vect,NA.char="-9", ploidy=1,ncode=**NCODE**)

library(ape)
library(genetics)
library(spdep)

Y = apply(W$tab,2,f1)

pca1 <- dudi.pca(Y, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 3) #calcul de la PCA

png("**current_dir**Eigenvalues.png",width = 1500, height = 1000, res=200)
barplot(pca1$eig[1:50], main = "Eigenvalues") #graphique des valeurs propres
dev.off()

png("**current_dir**PCA.png",width = 1500, height = 1000, res=200)
s.class(pca1$li, W$pop, sub = "PCA 1-2", csub = 2) #graphique de la PCA
add.scatter.eig(pca1$eig[1:5], nf = 3, xax = 1, yax = 2, posi = "top") #ajout des valeurs propres en cartouche sur le graphique précédent
dev.off()
"""

DAPCchange="""
#Kmeans
fcl.BIC <- find.clusters(W, n.pca=**PCARETAIN**,n.clust=**pop**, stat="BIC", n.iter=5000, n.start=30, scale=FALSE)
dapc <- dapc(W, pop=fcl.BIC$grp, n.pca=**PCARETAIN**,n.da=**PCARETAIN**, scale=FALSE, pca.select="nbEig")

#graphique de la DAPC
png("**current_dir**DAPC_K**pop**.png",width = 1500, height = 1000, res=200)
scatter(dapc, xax=1, yax=2,col=rainbow(length(levels(dapc$grp))))#, posi="topleft", bg="white",ratio=0.3, csub=1.2)
dev.off()

#assignation des individus aux differents groupes
png("**current_dir**K**pop**_groupe.png",width = 1500, height = 1000, res=200)
assignplot(dapc, only.grp=NULL, subset=NULL, cex.lab=.4, pch=3)
dev.off()

#recuperer les assignations des individus dans les differents clusters DAPC
write.csv2(dapc$posterior, file="**current_dir**K**pop**.csv")


png("**current_dir**K**pop**.png",width = 30000, height = 1500, res=200)
compoplot(dapc, only.grp=NULL,# affiche que le groupe selectionné
          subset=NULL,      # pour affichier un sous ensemble des données
          new.pred=NULL,
          col=c("red","green","blue","yellow","orange","cyan","burlywood3","darkorchid3","deeppink2","gray35","chartreuse"),    #couleurs
          lab=name,
          cex.names = 0.6,
          legend=TRUE,
          txt.leg=NULL,
          ncol=**pop**,
          posi=list(x = 50, y = 1.1),  # position de la legende
          cleg=.8,
          bg=transp("white"))
dev.off()


"""


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
#	start=time.clock()
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='run_DAPC.py', description='''This Programme run DAPC for k i to m''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display run_DAPC version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-m', '--mat', metavar="<filename>",type=existant_file, required=True, dest = 'matriceParam', help = 'matrice file path')
	filesreq.add_argument('-o', '--order', metavar="<filename>",type=existant_file, required=True, dest = 'orderMatriceParam', help = 'file with re-order name of matrice')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-pca', '--pcanum', metavar="<int>", required=False, default = "NULL", dest = 'pcaParam', help = 'Number value of PCA retains (default = NULL)')
	files.add_argument('-pi', '--popi', metavar="<int>", type = int, default=2, required=False, dest = 'nbpopiParam', help = 'Number of pop Min (default = 2)')
	files.add_argument('-pm', '--popm', metavar="<int>", type = int, default=10, required=False, dest = 'nbpopmParam', help = 'Number of pop Max (default = 10)')
	#filesreq.add_argument('-ind', '--indiv', metavar="<int>", type = int, required=False, dest = 'nbindParam', help = 'Number of indiv in matrice')
	#filesreq.add_argument('-ma', '--marker', metavar="<int>", type = int, required=False, dest = 'nbmarkParam', help = 'Number of makers in matrice')
	#filesreq.add_argument('-n', '--ncode', metavar="<int>", type = int, required=False, dest = 'ncodeParam', help = 'ncode in Matrice file')
	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#         Welcome in run_DAPC (Version " + version + ")           #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	matriceParam = relativeToAbsolutePath(args.matriceParam)
	nbpopiParam = args.nbpopiParam
	nbpopmParam = args.nbpopmParam
	pcaParam = args.pcaParam
	#nbindParam = args.nbindParam
	#nbmarkParam = args.nbmarkParam+1
	#ncodeParam = args.ncodeParam
	orderMatriceParam = relativeToAbsolutePath(args.orderMatriceParam)

	# ajoute à la variable current_dir le chemin ou est executer le script
	current_dir = os.path.dirname(os.path.abspath(__file__))

	#check all path for matrice
	# Pour l'utilisation sous Linux:
	if "linux" in sys.platform:
		current_dir += "/"
		if "/" not in matriceParam:
			print("ERROR: not all path for matrice file !!!!")
			exit()
	# Pour l'utilisation sous windows:
	if "win" in sys.platform and "darwin" not in sys.platform:
		current_dir += "\\"
		if "\\" not in matriceParam :
			print("ERROR: not all path for matrice file !!!!")
			exit()

	basename = orderMatriceParam.split("/")[-1].split(".")[0]
	workingDir = "/".join(orderMatriceParam.split("/")[:-1])+"/"

	# Test si min < Max
	if nbpopiParam < 2:
		print("ERROR: min must be > 2 !!!!")
		exit()
	if nbpopiParam >= nbpopmParam:
		print("ERROR: min > max  !!!!!")
		exit()

	# création d'un répertoire de sortie

	outputDir = workingDir+basename+"/"

	if os.path.exists(outputDir):
		print("Warning , folder "+outputDir+" already exist !!!!")
		print("Do you want to remove all analysis? (y/n)\n")
		inp = None
		while inp == None or inp != "y" or inp != "n" or inp != "yes" or inp != "no":
			inp = input()
			if inp == "y":
				print("Remove directory "+outputDir)
				rmtree(outputDir)
				os.mkdir(outputDir)
				break
			if inp == "n":
				print("Programm exit\n")
				exit()
	else:
		os.mkdir(outputDir)


	# charge l'ordre a refaire
	orderList = loadInListCol(orderMatriceParam, 0)

	# copie de la matrice dans un dico
	dicoMatrice = loadInDictLine(matriceParam)

	# Comptage du nombre d'individus, de markers et ncode:
	nbindParam = len(dicoMatrice.keys())-1

	fileMat = open(matriceParam,"r")
	header = fileMat.readline()
	nbmarkParam = len(header.split("\t"))
	header = " \t"+"\t".join(header.split("\t")[1:])


	nbcode = fileMat.readline().split("\t")[1]
	while nbcode == "-9":
		nbcode = fileMat.readline().split("\t")[1]
	ncodeParam = len(nbcode)
	fileMat.close()

	# ouverture du nouveau fichier trier
	with open(outputDir+basename+"_Reorder.tab","w") as reorderMatriceFile:
		reorderMatriceFile.write(header)
		for ind in orderList:
			if ind not in dicoMatrice.keys():
				print("ERROR: The individu %s define in label file was not in the matrice file !!! Exit programme" % ind)
				exit()
			line = dicoMatrice[ind].split("\t")[0]+"\t"+"\t".join(dicoMatrice[ind].split("\t")[1:]).replace("999","-9")
			reorderMatriceFile.write(line)

	print("Nb individus: %i" % nbindParam)
	print("Nb markers: %i" % (int(nbmarkParam)-1))
	print("ncodeParam: %i" % ncodeParam)

	#ouverture du script R
	Rscript = open(outputDir+basename+"_R_DAPC.R","w")


	# Ajout du path du fichier matrice dans DAPCfix
	# modifie Script R pour adapter aux parametres rentrés
	dictToReplace = {
	"**MAKERS**"	:		str(nbmarkParam),
	"**NCODE**"	:			str(ncodeParam),
	"**INDIV**"	:			str(nbindParam),
	"**PATHTOFILE**":		str(reorderMatriceFile.name),
	"**current_dir**"	:	str(outputDir)
	}

	DAPCfixModif = replace_all(dictToReplace, DAPCfix)

	#print(DAPCfixModif)
	Rscript.write(DAPCfixModif)

	for pop in range(nbpopiParam,nbpopmParam+1):
		#print(pop)
		popstr=str(pop)
		DAPCchange2 = DAPCchange.replace("**pop**",popstr).replace("**current_dir**",str(outputDir)).replace("**PCARETAIN**",str(pcaParam))
		#print(DAPCchange2)
		Rscript.write(DAPCchange2)

	Rscript.close()






		## Display a summary of the execution
	print("\n\nExecution summary:")

	print("  - Outputting \n\
\t-"+outputDir+basename+"_R_DAPC.R script created" )

	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
