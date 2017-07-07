#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
## @package GUI_DAPC.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################

import sys, os
#current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
##sys.path.insert(1,current_dir+'../modules/')
#sys.path.insert(1,current_dir+'../modules/')
#from MODULES_SEB import relativeToAbsolutePath, directory, existant_file, loadInList

# Python modules
import argparse
import shutil, re
import subprocess
from shutil import rmtree

from PyQt4.QtCore import QObject, SIGNAL
from PyQt4.QtGui import QFileDialog, QApplication, QIcon, QPixmap, QPlainTextEdit
from PyQt4 import uic
import syntax

##################################################
## Variables Globales
version = '1.0'
VERSION_DATE = '31/08/2016'


##################################################
## PARAM DAPC R
installPackageR = """
#!/usr/bin/Rscript --vanilla
# -*- coding: utf-8 -*-
# @author Sébastien RAVEL
# Build with GUI_DAPC.py

# more info on:
# http://nas-bgpi.myds.me:30002/projects/tutoriel/wiki/Lancer_la_DAPC?parent=Wiki

#If parckage not install please use command line before in R prompt:
#install.packages(c("ade4","adegenet", "ggplot2","ape", "genetics", "spdep" ))

"""

DAPCfix = """
# Mean Fonction use to remove NA value
f1 <- function(vec) {
  m <- mean(vec, na.rm = TRUE)
  vec[is.na(vec)] <- m
  return(vec)
}

myCol = c("red","green","blue","purple","orange","cyan","burlywood3","darkorchid3","deeppink2","gray35","chartreuse")

graphique = **GRAPH**

library(ade4)
library(adegenet)
library(ggplot2)

# Load data into dataframe
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


####################################################
## manual select value for PCA DA value and view BIC graph
####################################################
# For the first time, if you keep PCA value to NULL, change n.clust=2 to n.clust=NULL, then run line
# this use intercative mode to help choose the PCA retains value with graphic.
# then produce BIC graphic to choose number of pop optimal.

fcl.BIC <- find.clusters(W, n.pca=NULL,n.clust=NULL, stat="BIC", n.iter=5000, n.start=30, scale=FALSE)
dapc <- dapc(W, pop=fcl.BIC$grp, n.pca=NULL,n.da=NULL, scale=FALSE, pca.select="nbEig")

# When you are choose the PCA and DA value, add to variable or go to GUI_DAPC and rebuild R script with value

# Change value for all K
PCARETAINValue <- **PCARETAIN**
DARETAINValue <- **DARETAIN**

"""

DAPCchange = """

####################################################
## TEST K = **pop** value
####################################################
#Kmeans
fcl.BIC <- find.clusters(W, n.pca=PCARETAINValue,n.clust=**pop**, stat="BIC", n.iter=5000, n.start=30, scale=FALSE)
dapc <- dapc(W, pop=fcl.BIC$grp, n.pca=PCARETAINValue,n.da=DARETAINValue, scale=FALSE, pca.select="nbEig")

#graphique de la DAPC
png("**current_dir**DAPC_K**pop**.png",width = 1500, height = 1000, res=200)
if (graphique == 1)
{
	scatter(dapc, xax=1, yax=2,col=rainbow(length(levels(dapc$grp))))#, posi="topleft", bg="white",ratio=0.3, csub=1.2)
}
if (graphique == 2)
{
	scatter(dapc, bg="white", pch=17:27, cstar=0, col=myCol, scree.da=FALSE,scree.pca=FALSE)
}
if (graphique == 3)
{
	scatter(dapc,scree.pca=FALSE, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=myCol, solid=.4, cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:**pop**))
}
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
          col=c("red","green","blue","purple","orange","cyan","burlywood3","darkorchid3","deeppink2","gray35","chartreuse"),    #couleurs
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
## Functions

def resource_path(relative):
	if hasattr(sys, "_MEIPASS"):
		return os.path.join(sys._MEIPASS, relative)
	return os.path.join(relative)

def loadInListCol(filename, col):
	"""
	Load file in list() and then remove \\n at end of line

	:param filename: a file
	:type filename: file
	:rtype: list()
	:return: - list of row's file without \\n
	:warn: Use this function with small file !!! except more RAM are use and crash systeme.

	Example:
		>>> rows = loadInListCol(filename, 0)
		>>> rows
		["i like pears, but apples scare me","i like apples, but pears scare me","End of file"]
	"""

	liste = open(filename,"r").readlines()
	listgood=[line.rstrip().split("\t")[col] for line in liste]
	return listgood

def loadInDictLine(filename):
	"""
	Load file in Dict() and then remove \\n at end of line, then add first column in key of dict and valueare other column.

	:param filename: a file
	:type filename: file
	:rtype: dict()
	:return: - dict of row's file without \\n with key is first column and value list of other column
	:warn: Use this function with small file !!! except more RAM are use and crash systeme.

	Example:
		>>> dico = loadInDictLine(filename)
		>>> dico
		{
		"col1",[line1],
		"indiv1",[line2],
		"indiv2",[line3]
		}
	"""

	dicoOut={}
	with open(filename) as filein:
		for line in filein:
			tabLine = line.rstrip().split("\t")
			#print(tabLine[0], tabLine[1])
			if tabLine[0] not in dicoOut.keys():
				dicoOut[tabLine[0]] = line
	return dicoOut

def replace_all(repls, str):
	"""
	Function that take a dictionnary and text variable and return text variable with replace 'Key' from dictionnary with 'Value'.

	:param repls: a python dictionary
	:type repls: dict()
	:param str: a string where remplace some words
	:type str: str()
	:rtype: str()
	:return: - txt with replace 'Key' of dictionnary with 'Value' in the input txt

	Example:
		>>> text =  "i like apples, but pears scare me"
		>>> print(replace_all({"apple": "pear", "pear": "apple"}, text))
		i like pears, but apples scare me
	"""

	return re.sub('|'.join(re.escape(key) for key in repls.keys()),lambda k: repls[k.group(0)], str)

def existant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open by default.

	:param x: a file path
	:type x: str()
	:rtype: string
	:return: string

	"""
	if not os.path.exists(x):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))

	return x

def relativeToAbsolutePath(relative):
	"""	Return the absolutPath

	:param relative: a string path
	:type relative: string
	:rtype: string()
	:return: absolutePath
	:warn: need subprocess::check_output

	Example:
	>>>print(relative)
	../test
	>>> pathDirectory = relativeToAbsolutePath(relative)
	>>>print(pathDirectory)
	/home/sebastien/test

	"""
	from subprocess import check_output
	if relative[0] != "/":			# The relative path is a relative path, ie do not starts with /
		command = "readlink -m "+relative
		absolutePath = subprocess.check_output(command, shell=True).encode("utf-8").decode("utf-8").rstrip()
		return absolutePath
	else:						# Relative is in fact an absolute path, send a warning
		absolutePath = relative;
		return absolutePath


filename = './includes/dapc.ui'
myfontfile = resource_path(os.path.join(filename))
formDAPC,baseDAPC = uic.loadUiType(myfontfile)

def resource_path(relative_path):
	if hasattr(sys, '_MEIPASS'):
		return os.path.join(sys._MEIPASS, relative_path)
	return os.path.join(os.path.abspath("."), relative_path)

# ************************************* CLASSE DAPC Gestion de l'affichage et de la récupération de donnée ****************************************************
## @class DAPC
# @brief Classe principale, fenêtre principale
class DAPC( formDAPC, baseDAPC ):
	""" Classe principale qui est aussi la fenêtre principale de la GUI
	"""
	def __init__(self,app,parent=None):
		super(DAPC,self).__init__(parent)
		self.app = app
		self.initialiseVariables()
		self.createWidgets()


	def initialiseVariables(self):
		"""Initialise les variables par défaults"""

		self.PCAvalue = "NULL"
		self.DAvalue = "NULL"
		self.popMinValue = "2"
		self.popMaxValue = "10"
		self.matricePathFile = ""
		self.orderPathFile = ""
		self.workingDir = ""
		self.basename = ""
		self.rmOld = "False"
		self.DAPCfix = DAPCfix
		self.DAPCchange = DAPCchange
		self.expertMode = "False"
		self.graphType = "1"



		self.pathFileOut = ""

	def createWidgets(self):
		"""Mise en place du masquage des frames et des connections de l'interface"""
		self.ui = self
		self.ui.setupUi(self)
		self.ui.setFocus(True)

		## Initialisation des frames et bouttons
		self.ui.runPushButton.setDisabled(True)
		self.ui.expertFrame.hide()
		self.ui.frameRun.hide()
		#self.ui.progressbar.hide()
		self.setWindowIcon(QIcon(resource_path("./includes/icon.ico")))
		self.ui.statusbar.setStyleSheet("color: rgb(255, 107, 8);font: 14pt 'Arial Black';")

		pic = QPixmap(resource_path("./includes/graph"+self.graphType+".png"))
		pic = pic.scaled(298,298)
		self.ui.imgLabel.setPixmap(pic)

		self.resizeWindows()

		## Edition des connect:
		QObject.connect(self.ui.loadMatriceFilePushButton,SIGNAL("clicked()"),self.loadMatriceFile)
		QObject.connect(self.ui.loadOrderFilePushButton,SIGNAL("clicked()"),self.loadOrderFile)

		QObject.connect(self.ui.resetPushButton,SIGNAL("clicked()"),self.resetPress)
		QObject.connect(self.ui.runPushButton,SIGNAL("clicked()"),self.run)

		QObject.connect(self.ui.updateColorPushButton,SIGNAL("clicked()"),self.actualizeDAPCColor)

		QObject.connect(self.ui.rmOldCheckBox,SIGNAL("stateChanged(int)"),self.actualizeRmOld)
		QObject.connect(self.ui.expertCheckBox,SIGNAL("stateChanged(int)"),self.actualizeExpertMode)

		QObject.connect(self.ui.graphTypeComboBox,SIGNAL("currentIndexChanged(int)"),self.actualizeGraph)

		QObject.connect(self.ui.PCAlineEdit,SIGNAL("editingFinished()"),self.actualizePCA)
		QObject.connect(self.ui.DAlineEdit,SIGNAL("editingFinished()"),self.actualizeDA)
		QObject.connect(self.ui.popMinLineEdit,SIGNAL("editingFinished()"),self.actualizePopMin)
		QObject.connect(self.ui.popMaxLineEdit,SIGNAL("editingFinished()"),self.actualizePopMax)
		#QObject.connect(self.ui.DAPCfixPlainTextEdit,SIGNAL("blockCountChanged()"),self.actualizeDAPCFix)
		#QObject.connect(self.ui.DAPCchangePlainTextEdit,SIGNAL("blockCountChanged()"),self.actualizeDAPCChange)

	def resizeWindows(self):
		"""change la taille de fenetre"""
		# resize
		size = self.ui.centralwidget.sizeHint()
		self.ui.setGeometry(300,100,size.width(), size.height())


	def actualizeGraph(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.graphTypeComboBox.currentText())
		self.graphType = str(currentValue)
		pic = QPixmap(resource_path("./includes/graph"+self.graphType+".png"))
		pic = pic.scaled(298,298)
		self.ui.imgLabel.setPixmap(pic)


		#elif self.graphType in [2]:
			#self.ui.nuclOutLineEdit.setDisabled(True)

		#else:
			#self.ui.nuclOutLineEdit.setEnabled(True)

	def actualizeDAPCFix(self):
		"""change la valeur du choix quand changer"""
		#print("changeFix")
		currentValue = str(self.ui.DAPCfixPlainTextEdit.toPlainText().toUtf8())
		highlight = syntax.PythonHighlighter(self.ui.DAPCfixPlainTextEdit.document())
		self.ui.DAPCfixPlainTextEdit.setPlainText(currentValue.decode("utf-8"))

	def actualizeDAPCChange(self):
		"""change la valeur du choix quand changer"""
		#print("changeChange")
		currentValue = str(self.ui.DAPCchangePlainTextEdit.toPlainText().toUtf8())
		highlight = syntax.PythonHighlighter(self.ui.DAPCchangePlainTextEdit.document())
		self.ui.DAPCchangePlainTextEdit.setPlainText(currentValue.decode("utf-8"))

	def actualizeDAPCColor(self):
		"""change la valeur du choix quand changer"""
		self.actualizeDAPCFix()
		self.actualizeDAPCChange()



	def actualizeExpertMode(self):
		"""change la valeur du choix quand changer"""
		if self.ui.expertCheckBox.isChecked():
			self.expertMode = "True"
			self.ui.expertFrame.show()
			self.ui.DAPCfixPlainTextEdit.clear()
			self.ui.DAPCchangePlainTextEdit.clear()

			highlight = syntax.PythonHighlighter(self.ui.DAPCfixPlainTextEdit.document())
			highlight2 = syntax.PythonHighlighter(self.ui.DAPCchangePlainTextEdit.document())

			self.ui.DAPCfixPlainTextEdit.setPlainText(DAPCfix.decode("utf-8"))
			self.ui.DAPCchangePlainTextEdit.setPlainText(DAPCchange.decode("utf-8"))
		else:
			self.expertMode = "False"
			self.ui.expertFrame.hide()

	def actualizeRmOld(self):
		"""change la valeur du choix quand changer"""
		if self.ui.rmOldCheckBox.isChecked():
			self.rmOld = "True"
		else:
			self.rmOld = "False"

	def actualizePCA(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.PCAlineEdit.text())
		if currentValue != "" and currentValue != "NULL":
			try:
				int(currentValue)
				self.PCAvalue = currentValue
			except Exception as e:
				self.displayError(error = str(e)+"\n Please enter an integer value !!!!")
				self.PCAvalue = "NULL"
				self.ui.PCAlineEdit.setText("NULL")
		else:
			self.PCAvalue = "NULL"
			self.ui.PCAlineEdit.setText("NULL")

	def actualizeDA(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.DAlineEdit.text())
		if currentValue != "" and currentValue != "NULL":
			try:
				int(currentValue)
				self.DAvalue = currentValue
			except Exception as e:
				self.displayError(error = str(e)+"\n Please enter an integer value !!!!")
				self.DAvalue = "NULL"
				self.ui.DAlineEdit.setText("NULL")
		else:
			self.DAvalue = "NULL"
			self.ui.DAlineEdit.setText("NULL")

	def actualizePopMin(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.popMinLineEdit.text())
		try:
			if currentValue != "" and int(currentValue) >= 2 and int(currentValue) <= int(self.popMaxValue) :
				self.popMinValue = currentValue
			else:
				self.popMinValue = "2"
				self.ui.popMinLineEdit.setText("2")
		except Exception as e:
			self.displayError(error = str(e)+"\n Please enter an integer value !!!!")
			self.popMinValue = "2"
			self.ui.popMinLineEdit.setText("2")

	def actualizePopMax(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.popMaxLineEdit.text())
		try:
			if currentValue != "" and int(currentValue) >= int(self.popMinValue):
				self.popMaxValue = currentValue
			else:
				self.popMaxValue = "10"
				self.ui.popMaxLineEdit.setText("10")
		except Exception as e:
			self.displayError(error = str(e)+"\n Please enter an integer value !!!!")
			self.popMaxValue = "10"
			self.ui.popMaxLineEdit.setText("10")

	def loadMatriceFile(self):
		"""Méthode qui permet de charger un fichier et afficher dans le plainText"""
		filename = QFileDialog.getOpenFileName(self, caption="Load the File", directory=os.getcwd(), filter="Text files (*.txt *.tab);;All (*.*)")
		if filename == "" and self.matricePathFile == "":
			self.displayError(error = "\"%s\" is not a valid file \n" % (filename))
		else:
			self.matricePathFile = filename
			self.ui.matricePlainTextEdit.setPlainText(filename)
			self.actualizeRunBouton()

	def loadOrderFile(self):
		"""Méthode qui permet de charger un fichier et afficher dans le plainText"""
		filename = QFileDialog.getOpenFileName(self, caption="Load the File", directory=os.getcwd(), filter="Text files (*.txt *.tab);;All (*.*)")
		if filename == "" and self.orderPathFile == "":
			self.displayError(error = "\"%s\" is not a valid file \n" % (filename))
		else:
			self.orderPathFile = filename
			self.workingDir = "/".join(str(filename).split("/")[:-1])+"/"
			self.basename = filename.split("/")[-1].split(".")[0]
			self.pathFileOut = self.workingDir+self.basename+"/"
			if os.path.isdir(self.pathFileOut):
				self.displayError(error = "Warnnig folder %s already exist,\nPlease remove/rename before run new analysis or use checkbox\n" % (self.pathFileOut))
			self.ui.orderPlainTextEdit.setPlainText(filename)
			self.actualizeRunBouton()

	def actualizeRunBouton(self):
		"""de grise le bouton run si rempli matrice et order"""
		if self.matricePathFile != "" and self.orderPathFile != "":
			self.ui.runPushButton.setEnabled(True)
		else:
			self.ui.runPushButton.setDisabled(True)

	def resetPress(self):
		#"""reset PATH files"""
		self.ui.matricePlainTextEdit.clear()
		self.ui.orderPlainTextEdit.clear()
		self.ui.frameRun.hide()
		self.ui.runningPlainTextEdit.clear()

		self.ui.PCAlineEdit.setText("NULL")
		self.ui.DAlineEdit.setText("NULL")
		self.ui.popMinLineEdit.setText("2")
		self.ui.popMaxLineEdit.setText("10")

		self.ui.rmOldCheckBox.setChecked(False)
		self.ui.expertCheckBox.setChecked(False)

		#self.ui.progressbar.hide()
		self.ui.expertFrame.hide()
		self.ui.runPushButton.setDisabled(True)
		self.ui.loadMatriceFilePushButton.setEnabled(True)
		self.ui.loadOrderFilePushButton.setEnabled(True)
		self.ui.PCAlineEdit.setEnabled(True)
		self.ui.DAlineEdit.setEnabled(True)
		self.ui.popMinLineEdit.setEnabled(True)
		self.ui.popMaxLineEdit.setEnabled(True)
		self.ui.rmOldCheckBox.setEnabled(True)
		self.ui.expertCheckBox.setEnabled(True)
		self.ui.expertFrame.setEnabled(True)
		self.ui.graphTypeComboBox.setEnabled(True)


		self.graphTypeComboBox.setCurrentIndex(0)
		self.initialiseVariables()


	def run(self):
		try:
			warning = ""	# initialise le nombre d'erreur
			val = 0			# initialise le nombre d'erreur
			txtInfo = ""

			##grise les boutons pour pas relancer job
			self.ui.frameRun.show()
			self.ui.runPushButton.setDisabled(True)
			self.ui.loadMatriceFilePushButton.setDisabled(True)
			self.ui.loadOrderFilePushButton.setDisabled(True)
			self.ui.PCAlineEdit.setDisabled(True)
			self.ui.DAlineEdit.setDisabled(True)
			self.ui.popMinLineEdit.setDisabled(True)
			self.ui.popMaxLineEdit.setDisabled(True)
			self.ui.rmOldCheckBox.setDisabled(True)
			self.ui.expertCheckBox.setDisabled(True)
			self.ui.expertFrame.setDisabled(True)
			self.ui.graphTypeComboBox.setDisabled(True)
			if self.expertMode == "True":
				self.DAPCfix = str(self.ui.DAPCfixPlainTextEdit.toPlainText().toUtf8())
				self.DAPCchange = str(self.ui.DAPCchangePlainTextEdit.toPlainText().toUtf8())

			"""to run programme"""
			# création du pathout
			if os.path.isdir(self.pathFileOut):
				if self.rmOld == "True":
					shutil.rmtree(str(self.pathFileOut))
					os.mkdir(self.pathFileOut)
				else:
					warning += "Warnnig folder "+self.pathFileOut+" already exist,\nPlease remove/rename before run new analysis or use checkbox"
					raise Exception(warning)
			else:
				os.mkdir(self.pathFileOut)

			###############################################
			# code commun mode graphique ou interface
			###############################################

			# charge l'ordre a refaire
			self.orderList = loadInListCol(self.orderPathFile, 0)

			# copie de la matrice dans un dico
			self.dicoMatrice = loadInDictLine(self.matricePathFile)

			# Comptage du nombre d'individus, de markers et ncode:
			self.nbindParam = len(self.dicoMatrice.keys())-1

			if self.nbindParam != len(self.orderList):
				txtInfo += "WARNING: More individu in Matrice file (%s) than Order label file (%s)!!!\n" % (self.nbindParam, len(self.orderList))

			fileMat = open(self.matricePathFile,"r")
			header = fileMat.readline()
			self.nbmarkParam = len(header.split("\t"))
			header = " \t"+"\t".join(header.split("\t")[1:])


			nbcode = fileMat.readline().split("\t")[1]
			while nbcode == "-9":
				nbcode = fileMat.readline().split("\t")[1]
			self.ncodeParam = len(nbcode)
			fileMat.close()

			# ouverture du nouveau fichier trier
			with open(self.pathFileOut+self.basename+"_Reorder.tab","w") as reorderMatriceFile:
				reorderMatriceFile.write(header)
				for ind in self.orderList:
					if ind not in self.dicoMatrice.keys():
						error = "ERROR: The individu %s define in label file was not in the matrice file !!! Exit programme" % ind
						raise Exception(error)

					line = self.dicoMatrice[ind].split("\t")[0]+"\t"+"\t".join(self.dicoMatrice[ind].split("\t")[1:]).replace("999","-9")
					reorderMatriceFile.write(line)

			txtInfo += "Nb individus: %i\tNb markers: %i\tncodeParam: %i\tGraph type: %s\n" % (self.nbindParam,int(self.nbmarkParam)-1,self.ncodeParam, self.graphType)
			if args.cmdMode:
				pass
			else:
				self.ui.runningPlainTextEdit.setPlainText(txtInfo)

			#ouverture du script R
			Rscript = open(self.pathFileOut+self.basename+"_R_DAPC.R","w")
			Rscript.write(installPackageR)

			# Ajout du path du fichier matrice dans DAPCfix
			# modifie Script R pour adapter aux parametres rentrés
			dictToReplace = {
			"**MAKERS**"	:		str(self.nbmarkParam),
			"**NCODE**"	:			str(self.ncodeParam),
			"**INDIV**"	:			str(self.nbindParam),
			"**PATHTOFILE**":		str(reorderMatriceFile.name),
			"**current_dir**"	:	str(self.pathFileOut),
			"**GRAPH**"	:	str(self.graphType)
			}

			DAPCfixModif = replace_all(dictToReplace, self.DAPCfix)

			#print(DAPCfixModif)

			Rscript.write(DAPCfixModif)

			for pop in range(int(self.popMinValue),int(self.popMaxValue)+1):
				#print(pop)
				popstr=str(pop)
				DAPCchange2 = self.DAPCchange.replace("**pop**",popstr).replace("**current_dir**",str(self.pathFileOut)).replace("**PCARETAIN**",str(self.PCAvalue)).replace("**DARETAIN**",str(self.DAvalue))
				#print(DAPCchange2)
				Rscript.write(DAPCchange2)

			Rscript.close()
			self.ui.statusbar.showMessage(str("FINISH, script product on : %s" % self.pathFileOut),9600)

			txtInfo += "FINISH, script product on :\n %s" % (self.pathFileOut)

			if args.cmdMode:
				print(txtInfo)
			else:
				self.ui.runningPlainTextEdit.setPlainText(txtInfo)

			## si des erreurs:

		except Exception as e:
			self.displayError(error = e)



	def displayError(self, error):
		""" affiche les erreurs dans la zone de text"""
		if args.cmdMode:
				print(error)
		else:
			# Grise les cases pour ne pas relancer dessus et faire un reset
			self.ui.runPushButton.setDisabled(True)
			txtError = str(error)
			self.ui.statusbar.showMessage(txtError,7200)





def main():

	#print sys.argv+["-cmd", "-gnome-terminal"]
	nargv = sys.argv

	# instanciation des objets principaux
	app = QApplication(nargv)
	myapp = DAPC(app)

	myapp.showMinimized()
	# les .app sous macos nécessitent cela pour que l'appli s'affiche en FG
	if "darwin" in sys.platform and ".app/" in sys.argv[0]:
		myapp.raise_()

	# lancement de la boucle Qt
	sys.exit(app.exec_())


def cmd():
	#print sys.argv+["-cmd", "-gnome-terminal"]
	nargv = sys.argv
	# instanciation des objets principaux
	app = QApplication(nargv)
	myapp = DAPC(app)

	# Load info arguments to Class
	# matrice file
	myapp.matricePathFile = relativeToAbsolutePath(args.matriceParam)
	# order file
	myapp.orderPathFile = relativeToAbsolutePath(args.orderMatriceParam)
	# PCA value
	myapp.PCAvalue = args.pcaParam
	# DA value
	myapp.DAvalue = args.daParam
	# pop min value
	myapp.popMinValue = args.nbpopiParam
	# pop max value
	myapp.popMaxValue = args.nbpopmParam
	# pgraph value
	myapp.graphType = args.graphParam

	# working dir path
	workingDir = "/".join(relativeToAbsolutePath(args.orderMatriceParam).encode("utf-8").split("/")[:-1])+"/"
	myapp.workingDir = workingDir.encode("utf-8")
	# basename
	basename = relativeToAbsolutePath(args.orderMatriceParam).encode("utf-8").split("/")[-1].split(".")[0]
	myapp.basename = basename.encode("utf-8")
	# pathFileOut dir path
	pathFileOut = workingDir+basename+"/"
	myapp.pathFileOut = pathFileOut.encode("utf-8")

	# rm old folder
	myapp.rmOld = args.rmOldParam


	# Run programme
	myapp.run()


if __name__ == '__main__':

	# Parameters recovery
	parser = argparse.ArgumentParser(prog='GUI_DAPC.py', description='''This Programme open GUI to produce DAPC script.\n
																				#If use on cluster you can run in commande line with option -c and args''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display GUI_DAPC.py version number and exit')
	parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running if -c use')
	filesReq.add_argument('-c', '--cmd', action='store_true', dest = 'cmdMode', help = 'If used, programme run in CMD without interface')
	filesReq.add_argument('-m', '--mat', metavar="<filename>",type=existant_file, required=False, dest = 'matriceParam', help = 'matrice file path')
	filesReq.add_argument('-o', '--order', metavar="<filename>",type=existant_file, required=False, dest = 'orderMatriceParam', help = 'file with re-order name of matrice')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-pca', '--pcanum', metavar="<int>", required=False, default = "NULL", dest = 'pcaParam', help = 'Number value of PCA retains (default = NULL)')
	files.add_argument('-da', '--danum', metavar="<int>", required=False, default = "NULL", dest = 'daParam', help = 'Number value of DA retains (default = NULL)')
	files.add_argument('-pi', '--popi', metavar="<int>", type = int, default=2, required=False, dest = 'nbpopiParam', help = 'Number of pop Min (default = 2)')
	files.add_argument('-pm', '--popm', metavar="<int>", type = int, default=10, required=False, dest = 'nbpopmParam', help = 'Number of pop Max (default = 10)')	## Check parameters
	files.add_argument('-r', '--rm', metavar="<True/False>", choices=("True","False"), default="False", required=False, dest = 'rmOldParam', help = 'if directory exist remove (default = False)')	## Check parameters
	files.add_argument('-g', '--graph', metavar="<1/2/3>", choices=("1","2","3"), default="1", required=False, dest = 'graphParam', help = 'type of graph (default = 1)')	## Check parameters
	args = parser.parse_args()

	if args.cmdMode:
		if args.matriceParam == "" or args.orderMatriceParam == "":
			print "ERROR: You must enter require arguments"
			exit()
		cmd()
	else:
		# run interface
		main()
