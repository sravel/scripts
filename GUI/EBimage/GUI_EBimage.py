#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
## @package GUI_EBimage.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################

import sys, os

# Python modules
import argparse
import shutil, re
import subprocess
from shutil import rmtree

from PyQt5.QtCore import QObject, SIGNAL
from PyQt5.QtGui import QFileDialog, QApplication, QIcon, QPixmap, QPlainTextEdit, QMessageBox
from PyQt5 import uic
import glob

import syntax

from rpy2.robjects.packages import importr
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

#exit()
##################################################
## Variables Globales
version = '1.0'
VERSION_DATE = '02/10/2017'


##################################################
## PARAM EBimage R
installPackageR = """
#!/usr/bin/Rscript --vanilla
# -*- coding: utf-8 -*-
# @author Sébastien RAVEL
# Build with GUI_EBimage.py

# more info on:
# http://nas-bgpi.myds.me:30002/projects/tutoriel/wiki/Lancer_la_EBimage?parent=Wiki

#If parckage not install please use command line before in R prompt:
#install.packages(c("ade4","adegenet", "ggplot2","ape", "genetics", "spdep" ))

"""

EBimageCallibration = """

"""

EBimageAnalysis = """

"""


##################################################
## Functions

def resource_path(relative):
	if hasattr(sys, "_MEIPASS"):
		return os.path.join(sys._MEIPASS, relative)
	return os.path.join(relative)


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

filename = './includes/EBimage.ui'
myfontfile = resource_path(os.path.join(filename))
formEBimage,baseEBimage = uic.loadUiType(myfontfile)

def resource_path(relative_path):
	if hasattr(sys, '_MEIPASS'):
		return os.path.join(sys._MEIPASS, relative_path)
	return os.path.join(os.path.abspath("."), relative_path)

# ************************************* CLASSE EBimage Gestion de l'affichage et de la récupération de donnée ****************************************************
## @class EBimage
# @brief Classe principale, fenêtre principale
class EBimage( formEBimage, baseEBimage ):
	""" Classe principale qui est aussi la fenêtre principale de la GUI
	"""
	def __init__(self,app,parent=None):
		super(EBimage,self).__init__(parent)
		self.app = app
		self.ui = self
		self.ui.setupUi(self)
		self.ui.setFocus(True)
		self.initialiseVariables()
		self.createWidgets()


	def initialiseVariables(self):
		"""Initialise les variables par défaults"""

		# callibration
		self.dicoFoldersCallibration =	{"leaf": None,
							 "symptom": None,
							 "background": None,
							}
							#"cat1": None,
							#"cat2": None,

		self.futherCat = 0
		self.dicoObjectOpenLineEditCallibration = {
						"leaf": self.ui.leafOpenLineEdit,
						"symptom": self.ui.symptomOpenLineEdit,
						"background": self.ui.backgroundOpenLineEdit,
					 }
						#"cat1": self.ui.category1OpenLineEdit,
						#"cat2": self.ui.category2OpenLineEdit
		self.dicoObjectOpenPushButtonCallibration = {
						"leaf": self.ui.loadLeafPushButton,
						"symptom": self.ui.loadSymptomPushButton,
						"background": self.ui.loadBackgroundPushButton,
					 }
						#"cat1": self.ui.loadCategory1PushButton,
						#"cat2": self.ui.loadCategoryPushButton


		# Analysis
		self.leafMinSize = 20
		self.symptomMinSize = 0
		self.lesionBorder = ""
		self.inputFolder = ""
		self.outputFolder = ""
		#self.rmOld = "False"

		self.EBimageCallibration = EBimageCallibration
		self.EBimageAnalysis = EBimageAnalysis

	def createWidgets(self):
		"""Mise en place du masquage des frames et des connections de l'interface"""
		## Initialisation des frames et bouttons
		self.ui.runAnalysisPushButton.setDisabled(True)
		self.ui.runCallibrationPushButton.setDisabled(True)
		self.ui.runCallibrationFrame.hide()
		self.ui.futherFrame.hide()
		self.ui.cat1Frame.hide()
		self.ui.cat2Frame.hide()

		self.setWindowIcon(QIcon(resource_path("./includes/icon.ico")))
		self.ui.statusbar.setStyleSheet("color: rgb(255, 107, 8);font: 12pt 'Arial Black';")

		#pic = QPixmap(resource_path("./includes/graph"+self.graphType+".png"))
		#pic = pic.scaled(298,298)
		#self.ui.imgLabel.setPixmap(pic)

		#self.resizeWindows()

		## Edition des connect:
		QObject.connect(self.ui.loadCategory1PushButton,SIGNAL("clicked()"),lambda: self.loadFolder(inputCat = "cat1", methode="clicked"))
		QObject.connect(self.ui.category1OpenLineEdit,SIGNAL("editingFinished()"),lambda: self.loadFolder(inputCat = "cat1", methode="write"))

		QObject.connect(self.ui.loadCategoryPushButton,SIGNAL("clicked()"),lambda: self.loadFolder(inputCat = "cat2", methode="clicked"))
		QObject.connect(self.ui.category2OpenLineEdit,SIGNAL("editingFinished()"),lambda: self.loadFolder(inputCat = "cat2", methode="write"))

		QObject.connect( self.ui.loadLeafPushButton,SIGNAL("clicked()"),lambda: self.loadFolder(inputCat = "leaf", methode="clicked"))
		QObject.connect(self.ui.leafOpenLineEdit,SIGNAL("editingFinished()"),lambda: self.loadFolder(inputCat = "leaf", methode="write"))

		QObject.connect(self.ui.loadBackgroundPushButton,SIGNAL("clicked()"),lambda: self.loadFolder(inputCat = "background", methode="clicked"))
		QObject.connect(self.ui.backgroundOpenLineEdit,SIGNAL("editingFinished()"),lambda: self.loadFolder(inputCat = "background", methode="write"))

		QObject.connect(self.ui.loadSymptomPushButton,SIGNAL("clicked()"),lambda: self.loadFolder(inputCat = "symptom", methode="clicked"))
		QObject.connect(self.ui.symptomOpenLineEdit,SIGNAL("editingFinished()"),lambda: self.loadFolder(inputCat = "symptom", methode="write"))

		QObject.connect(self.ui.futherCategoriesComboBox,SIGNAL("currentIndexChanged(int)"),self.actualizeFutherCat)

		QObject.connect(self.ui.runCallibrationPushButton,SIGNAL("clicked()"),self.run)


		#QObject.connect(self.ui.resetPushButton,SIGNAL("clicked()"),self.resetPress)

		#QObject.connect(self.ui.updateColorPushButton,SIGNAL("clicked()"),self.actualizeEBimageColor)

		#QObject.connect(self.ui.rmOldCheckBox,SIGNAL("stateChanged(int)"),self.actualizeRmOld)
		#QObject.connect(self.ui.expertCheckBox,SIGNAL("stateChanged(int)"),self.actualizeExpertMode)

		#QObject.connect(self.ui.graphTypeComboBox,SIGNAL("currentIndexChanged(int)"),self.actualizeGraph)

		#QObject.connect(self.ui.PCAlineEdit,SIGNAL("editingFinished()"),self.actualizePCA)
		#QObject.connect(self.ui.DAlineEdit,SIGNAL("editingFinished()"),self.actualizeDA)
		#QObject.connect(self.ui.popMinLineEdit,SIGNAL("editingFinished()"),self.actualizePopMin)
		#QObject.connect(self.ui.popMaxLineEdit,SIGNAL("editingFinished()"),self.actualizePopMax)
		#QObject.connect(self.ui.EBimagefixPlainTextEdit,SIGNAL("blockCountChanged()"),self.actualizeEBimageFix)
		#QObject.connect(self.ui.EBimagechangePlainTextEdit,SIGNAL("blockCountChanged()"),self.actualizeEBimageChange)

	def resizeWindows(self):
		"""change la taille de fenetre"""
		# resize
		size = self.ui.centralwidget.sizeHint()
		self.ui.setGeometry(800,600,size.width(), size.height())


	def loadFolder(self,inputCat = None, methode = None):
		"""Méthode qui permet de charger un fichier et afficher dans le plainText"""
		directoryToOpen = os.getcwd()
		for key, folder in self.dicoFoldersCallibration.items():
			if folder != None:
				print(key, folder)
				directoryToOpen = folder
				break

		if methode == "clicked":
			pathdir = QFileDialog.getExistingDirectory(self, caption="Load the directory "+inputCat, directory=directoryToOpen)
			self.dicoObjectOpenLineEditCallibration[inputCat].setText(pathdir)
		elif methode == "write":
			pathdir = str(self.dicoObjectOpenLineEditCallibration[inputCat].text())

		if pathdir != "" and os.path.isdir(pathdir):
			self.dicoFoldersCallibration[inputCat] = pathdir
		else:
			self.dicoFoldersCallibration[inputCat] = ""
			self.dicoObjectOpenLineEditCallibration[inputCat].setText("")
			self.displayError(error = "\"%s\" is not a valid Path \n" % (pathdir))
		self.actualizeRunBouton()
		self.autocomplete_Folder()
		#print(self.dicoFoldersCallibration[inputCat])

	def autocomplete_Folder(self):
		"""autocomplet les répertoires si un déja remplit"""
		directoryToOpen = None
		for key, folder in self.dicoFoldersCallibration.items():
			if folder != None:
				directoryToOpen = "/".join(str(folder).split("/")[:-1])
		if directoryToOpen != None:
			for root, dirs, files in os.walk(directoryToOpen):
				#print(root, dirs, files)
				for inputCat in self.dicoFoldersCallibration.keys():
					if inputCat in dirs:
						pathdir = "{}/{}".format(root,inputCat)
						self.dicoFoldersCallibration[inputCat] = pathdir
						self.dicoObjectOpenLineEditCallibration[inputCat].setText(pathdir)
		self.actualizeRunBouton()


	def actualizeFutherCat(self):
		"""masque les autres catégorie quand changer"""
		oldState = self.futherCat
		self.futherCat  = str(self.ui.futherCategoriesComboBox.currentText())
		if self.futherCat == "0":
			self.ui.cat1Frame.hide()
			self.ui.cat2Frame.hide()
			self.dicoFoldersCallibration["cat1"] = None
			self.dicoFoldersCallibration["cat2"] = None
			self.dicoObjectOpenLineEditCallibration["cat1"].setText("")
			self.dicoObjectOpenLineEditCallibration["cat2"].setText("")
			self.displayError(error = "Reset Futher categories path !! \n")
		elif self.futherCat == "1":
			self.ui.cat1Frame.show()
			self.ui.cat2Frame.hide()
			if oldState == "2":
				self.dicoObjectOpenLineEditCallibration["cat2"].setText("")
				self.displayError(error = "Reset categories 2 path !! \n")
			self.dicoFoldersCallibration["cat2"] = None
		elif self.futherCat == "2":
			self.ui.cat1Frame.show()
			self.ui.cat2Frame.show()
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
		if self.futherCat == 0:
			for key, folder in self.dicoFoldersCallibration.items():
				#print(key, folder)
				if folder != None :
					self.ui.runCallibrationPushButton.setEnabled(True)
				elif key not in ["cat1", "cat2"]:
					self.ui.runCallibrationPushButton.setDisabled(True)
					break
		else:
			if self.futherCat == "1" and self.dicoFoldersCallibration["cat1"] != None:
				self.ui.runCallibrationPushButton.setEnabled(True)
			elif self.futherCat == "2" and self.dicoFoldersCallibration["cat1"] != None and self.dicoFoldersCallibration["cat2"] != None:
				self.ui.runCallibrationPushButton.setEnabled(True)
			else:
				self.ui.runCallibrationPushButton.setDisabled(True)


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
		print("RUN")
		try:
			warning = ""	# initialise le nombre d'erreur
			val = 0			# initialise le nombre d'erreur
			txtInfo = ""

			# import R's "base" package
			base = importr('base')

			with open("fonctions_apprentissage.r", "r", encoding="utf-8") as apprentissageRopen:
				apprentissage = "".join(apprentissageRopen.readlines())
			print(apprentissage)

			apprentissage = SignatureTranslatedAnonymousPackage(apprentissage, "apprentissage")

			#path_sample = "/media/sebastien/Bayer/ScriptsSEB/scripts/GUI/EBimage/AnalyseImagesV2/Samples/5583"
			path_sample = "/".join(str(self.dicoFoldersCallibration["leaf"]).split("/")[:-1])
			print(path_sample)
			print(dir(apprentissage))
			result , pathRdataFile= apprentissage.apprentissage(path_sample).r_repr().replace('"','').replace("c(","").replace(")","").split(",")
			print(result, pathRdataFile)
			if result == "1":
				reply = QMessageBox.question(parent=self, title='Attention', text='File will be overwritten.\nDo you still want to proceed?', buttons=QMessageBox.Yes | QMessageBox.No, defaultButton=QMessageBox.No)
				if reply == QMessageBox.Yes:
					print("OK")
				self.callibrationFileOpenLineEdit.setText(pathRdataFile)
			else:

				print("BAD")


			##grise les boutons pour pas relancer job
			#self.ui.frameRun.show()
			#self.ui.runPushButton.setDisabled(True)
			#self.ui.loadMatriceFilePushButton.setDisabled(True)
			#self.ui.loadOrderFilePushButton.setDisabled(True)
			#self.ui.PCAlineEdit.setDisabled(True)
			#self.ui.DAlineEdit.setDisabled(True)
			#self.ui.popMinLineEdit.setDisabled(True)
			#self.ui.popMaxLineEdit.setDisabled(True)
			#self.ui.rmOldCheckBox.setDisabled(True)
			#self.ui.expertCheckBox.setDisabled(True)
			#self.ui.expertFrame.setDisabled(True)
			#self.ui.graphTypeComboBox.setDisabled(True)
			#if self.expertMode == "True":
				#self.EBimagefix = str(self.ui.EBimagefixPlainTextEdit.toPlainText().toUtf8())
				#self.EBimagechange = str(self.ui.EBimagechangePlainTextEdit.toPlainText().toUtf8())

			#"""to run programme"""
			## création du pathout
			#if os.path.isdir(self.pathFileOut):
				#if self.rmOld == "True":
					#shutil.rmtree(str(self.pathFileOut))
					#os.mkdir(self.pathFileOut)
				#else:
					#warning += "Warnnig folder "+self.pathFileOut+" already exist,\nPlease remove/rename before run new analysis or use checkbox"
					#raise Exception(warning)
			#else:
				#os.mkdir(self.pathFileOut)

			################################################
			## code commun mode graphique ou interface
			################################################

			## charge l'ordre a refaire
			#self.orderList = loadInListCol(self.orderPathFile, 0)

			## copie de la matrice dans un dico
			#self.dicoMatrice = loadInDictLine(self.matricePathFile)

			## Comptage du nombre d'individus, de markers et ncode:
			#self.nbindParam = len(self.dicoMatrice.keys())-1

			#if self.nbindParam != len(self.orderList):
				#txtInfo += "WARNING: More individu in Matrice file (%s) than Order label file (%s)!!!\n" % (self.nbindParam, len(self.orderList))

			#fileMat = open(self.matricePathFile,"r")
			#header = fileMat.readline()
			#self.nbmarkParam = len(header.split("\t"))
			#header = " \t"+"\t".join(header.split("\t")[1:])


			#nbcode = fileMat.readline().split("\t")[1]
			#while nbcode == "-9":
				#nbcode = fileMat.readline().split("\t")[1]
			#self.ncodeParam = len(nbcode)
			#fileMat.close()

			## ouverture du nouveau fichier trier
			#with open(self.pathFileOut+self.basename+"_Reorder.tab","w") as reorderMatriceFile:
				#reorderMatriceFile.write(header)
				#for ind in self.orderList:
					#if ind not in self.dicoMatrice.keys():
						#error = "ERROR: The individu %s define in label file was not in the matrice file !!! Exit programme" % ind
						#raise Exception(error)

					#line = self.dicoMatrice[ind].split("\t")[0]+"\t"+"\t".join(self.dicoMatrice[ind].split("\t")[1:]).replace("999","-9")
					#reorderMatriceFile.write(line)

			#txtInfo += "Nb individus: %i\tNb markers: %i\tncodeParam: %i\tGraph type: %s\n" % (self.nbindParam,int(self.nbmarkParam)-1,self.ncodeParam, self.graphType)
			#if args.cmdMode:
				#pass
			#else:
				#self.ui.runningPlainTextEdit.setPlainText(txtInfo)

			##ouverture du script R
			#Rscript = open(self.pathFileOut+self.basename+"_R_EBimage.R","w")
			#Rscript.write(installPackageR)

			## Ajout du path du fichier matrice dans EBimagefix
			## modifie Script R pour adapter aux parametres rentrés
			#dictToReplace = {
			#"**MAKERS**"	:		str(self.nbmarkParam),
			#"**NCODE**"	:			str(self.ncodeParam),
			#"**INDIV**"	:			str(self.nbindParam),
			#"**PATHTOFILE**":		str(reorderMatriceFile.name),
			#"**current_dir**"	:	str(self.pathFileOut),
			#"**GRAPH**"	:	str(self.graphType)
			#}

			#EBimagefixModif = replace_all(dictToReplace, self.EBimagefix)

			##print(EBimagefixModif)

			#Rscript.write(EBimagefixModif)

			#for pop in range(int(self.popMinValue),int(self.popMaxValue)+1):
				##print(pop)
				#popstr=str(pop)
				#EBimagechange2 = self.EBimagechange.replace("**pop**",popstr).replace("**current_dir**",str(self.pathFileOut)).replace("**PCARETAIN**",str(self.PCAvalue)).replace("**DARETAIN**",str(self.DAvalue))
				##print(EBimagechange2)
				#Rscript.write(EBimagechange2)

			#Rscript.close()
			#self.ui.statusbar.showMessage(str("FINISH, script product on : %s" % self.pathFileOut),9600)

			#txtInfo += "FINISH, script product on :\n %s" % (self.pathFileOut)

			#if args.cmdMode:
				#print(txtInfo)
			#else:
				#self.ui.runningPlainTextEdit.setPlainText(txtInfo)

			## si des erreurs:

		except Exception as e:
			self.displayError(error = e)



	def displayError(self, error):
		""" affiche les erreurs dans la zone de text"""
		if args.cmdMode:
				print(error)
		else:
			# Grise les cases pour ne pas relancer dessus et faire un reset
			#self.ui.runPushButton.setDisabled(True)
			print(error)
			txtError = str(error)
			self.ui.statusbar.showMessage(txtError,7200)





def main():

	#print sys.argv+["-cmd", "-gnome-terminal"]
	nargv = sys.argv

	# instanciation des objets principaux
	app = QApplication(nargv)
	myapp = EBimage(app)

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
	myapp = EBimage(app)

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
	parser = argparse.ArgumentParser(prog='GUI_EBimage.py', description='''This Programme open GUI to produce EBimage script.\n
																				#If use on cluster you can run in commande line with option -c and args''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display GUI_EBimage.py version number and exit')
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
			print("ERROR: You must enter require arguments")
			exit()
		cmd()
	else:
		# run interface
		main()
