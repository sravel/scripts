#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
## @package GUI_structure.py
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
## PARAM Structure

qsubtxt = "qsub -V -q long.q -cwd "

mainparams = """\
// KEY PARAMETERS FOR THE PROGRAM structure.  YOU WILL NEED TO SET THESE
// IN ORDER TO RUN THE PROGRAM.  VARIOUS OPTIONS CAN BE ADJUSTED IN THE
// FILE extraparams.

//"(int)" means that this takes an integer value.
//"(B)"   means that this variable is Boolean
//        (ie insert 1 for True, and 0 for False)
//"(str)" means that this is a string (but not enclosed in quotes!)

// Basic Program Parameters

#define MAXPOPS   **POP**      // (int) number of populations assumed
#define BURNIN    100000       // (int) length of burnin period
#define NUMREPS   300000       // (int) number of MCMC reps after burnin

// Input/Output files

#define INFILE   **FILEIN**  // (str) name of input data file
#define OUTFILE  **FILEOUT**  //(str) name of output data file

//Data file format

#define NUMINDS    **nIndiv**    // (int) number of diploid individuals in data file
#define NUMLOCI    **nLoci**    // (int) number of loci in data file
#define PLOIDY       1    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line


#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   0     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says
                              whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier

#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data
                             before the genotype data start.

#define MARKERNAMES      1 // (B) data file contains row of marker names
#define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs)
                            // and a row to indicate which alleles are recessive
#define MAPDISTANCES     0  // (B) data file contains row of map distances
                            // between loci


// Advanced data file options

#define PHASED           0 // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO        0 // (B) the data for each individual contains a line
                                 // indicating phase (linkage model)
#define MARKOVPHASE      0 // (B) the phase info follows a Markov model.
#define NOTAMBIGUOUS  -9 // (int) for use in some analyses of polyploid data

"""

extraparams = """"""

##################################################
## Functions

def resource_path(relative):
	if hasattr(sys, "_MEIPASS"):
		return os.path.join(sys._MEIPASS, relative)
	return os.path.join(relative)

class printCol():
	"""
	printCol() CLASS
	================

	Classe qui ajoute des méthodes à print pour afficher de la couleur

	Example:

	>>> printCol.red("j'affiche en rouge")
	j'affiche en rouge

	"""

	RED = '\033[91m'
	GREEN = '\033[92m'
	YELLOW = '\033[93m'
	LIGHT_PURPLE = '\033[94m'
	PURPLE = '\033[95m'
	END = '\033[0m'

	@classmethod
	def red(cls, s):
		print(cls.RED + s + cls.END)

	@classmethod
	def green(cls, s):
		print(cls.GREEN + s + cls.END)

	@classmethod
	def yellow(cls, s):
		print(cls.YELLOW + s + cls.END)

	@classmethod
	def lightPurple(cls, s):
		print(cls.LIGHT_PURPLE + s + cls.END)

	@classmethod
	def purple(cls, s):
		print(cls.PURPLE + s + cls.END)

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


filename = './includes/structure.ui'
myfontfile = resource_path(os.path.join(filename))
formSTRUCTURE,baseSTRUCTURE = uic.loadUiType(myfontfile)

def resource_path(relative_path):
	if hasattr(sys, '_MEIPASS'):
		return os.path.join(sys._MEIPASS, relative_path)
	return os.path.join(os.path.abspath("."), relative_path)

# ************************************* CLASSE STRUCTURE Gestion de l'affichage et de la récupération de donnée ****************************************************
## @class STRUCTURE
# @brief Classe principale, fenêtre principale
class STRUCTURE( formSTRUCTURE, baseSTRUCTURE ):
	""" Classe principale qui est aussi la fenêtre principale de la GUI
	"""
	def __init__(self,app,parent=None):
		super(STRUCTURE,self).__init__(parent)
		self.app = app
		self.initialiseVariables()
		self.createWidgets()


	def initialiseVariables(self):
		"""Initialise les variables par défaults"""

		# matrice
		self.matricePathFile = "NULL"
		# info tab
		self.nbIndiv = "NULL"
		self.nbMarker = "NULL"
		# K values
		self.popMinValue = "1"
		self.popMaxValue = "10"
		# rep values
		self.repMinValue = "1"
		self.repMaxValue = "10"

		# script build param
		self.workingDir = "NULL"
		self.pathFileOut = "NULL"
		self.basename = "NULL"
		self.rmOld = "False"
		self.mainparams = mainparams
		self.extraparams = extraparams
		self.outputPref = "NULL"


	def createWidgets(self):
		"""Mise en place du masquage des frames et des connections de l'interface"""
		self.ui = self
		self.ui.setupUi(self)
		self.ui.setFocus(True)

		## Initialisation des frames et bouttons
		self.ui.runPushButton.setDisabled(True)
		self.ui.frameRun.hide()

		# add to mainparams and extraparams
		self.actualizeSTRUCTUREColor()

		#self.ui.progressbar.hide()
		self.setWindowIcon(QIcon(resource_path("./includes/icon.ico")))
		self.ui.statusbar.setStyleSheet("color: rgb(255, 107, 8);font: 10pt 'Arial Black';")

		self.resizeWindows()

		## Edition des connect:
		QObject.connect(self.ui.loadMatriceFilePushButton,SIGNAL("clicked()"),self.loadMatriceFile)

		QObject.connect(self.ui.resetPushButton,SIGNAL("clicked()"),self.resetPress)
		QObject.connect(self.ui.runPushButton,SIGNAL("clicked()"),self.run)

		QObject.connect(self.ui.updateColorPushButton,SIGNAL("clicked()"),self.actualizeSTRUCTUREColor)

		QObject.connect(self.ui.rmOldCheckBox,SIGNAL("stateChanged(int)"),self.actualizeRmOld)

		QObject.connect(self.ui.indivLineEdit,SIGNAL("editingFinished()"),self.actualizeIndiv)
		QObject.connect(self.ui.markerLineEdit,SIGNAL("editingFinished()"),self.actualizeMarker)


		QObject.connect(self.ui.popMinLineEdit,SIGNAL("editingFinished()"),self.actualizePopMin)
		QObject.connect(self.ui.popMaxLineEdit,SIGNAL("editingFinished()"),self.actualizePopMax)

		QObject.connect(self.ui.repMinLineEdit,SIGNAL("editingFinished()"),self.actualizeRepMin)
		QObject.connect(self.ui.repMaxLineEdit,SIGNAL("editingFinished()"),self.actualizeRepMax)

		#QObject.connect(self.ui.mainparamsPlainTextEdit,SIGNAL("blockCountChanged()"),self.actualizeMainparams)
		#QObject.connect(self.ui.extraparamsPlainTextEdit,SIGNAL("blockCountChanged()"),self.actualizeExtraparams)

	def resizeWindows(self):
		"""change la taille de fenetre"""
		# resize
		size = self.ui.centralwidget.sizeHint()
		self.ui.setGeometry(300,200,size.width(), size.height())


	def actualizeMainparams(self):
		"""change la valeur du choix quand changer"""
		#print("changeFix")
		currentValue = str(self.ui.mainparamsPlainTextEdit.toPlainText().toUtf8())
		if currentValue == "":
			currentValue = mainparams
		highlight = syntax.PythonHighlighter(self.ui.mainparamsPlainTextEdit.document())
		self.ui.mainparamsPlainTextEdit.setPlainText(currentValue.decode("utf-8"))

	def actualizeExtraparams(self):
		"""change la valeur du choix quand changer"""
		#print("changeChange")
		currentValue = str(self.ui.extraparamsPlainTextEdit.toPlainText().toUtf8())
		if currentValue == "":
			currentValue = extraparams
		highlight = syntax.PythonHighlighter(self.ui.extraparamsPlainTextEdit.document())
		self.ui.extraparamsPlainTextEdit.setPlainText(currentValue.decode("utf-8"))

	def actualizeSTRUCTUREColor(self):
		"""change la valeur du choix quand changer"""
		self.actualizeMainparams()
		self.actualizeExtraparams()

	def actualizeRmOld(self):
		"""change la valeur du choix quand changer"""
		if self.ui.rmOldCheckBox.isChecked():
			self.rmOld = "True"
		else:
			self.rmOld = "False"

	def actualizeIndiv(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.indivLineEdit.text())
		if currentValue != "" and currentValue != "NULL":
			try:
				int(currentValue)
				self.nbIndiv = currentValue
				self.actualizeRunBouton()
			except Exception as e:
				self.displayError(error = str(e)+"\n Please enter an integer value !!!!")
				self.nbIndiv = "NULL"
				self.ui.indivLineEdit.setText("NULL")
		else:
			self.nbIndiv = "NULL"
			self.ui.indivLineEdit.setText("NULL")

	def actualizeMarker(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.markerLineEdit.text())
		if currentValue != "" and currentValue != "NULL":
			try:
				int(currentValue)
				self.nbMarker = currentValue
				self.actualizeRunBouton()
			except Exception as e:
				self.displayError(error = str(e)+"\n Please enter an integer value !!!!")
				self.nbMarker = "NULL"
				self.ui.markerLineEdit.setText("NULL")
		else:
			self.nbMarker = "NULL"
			self.ui.markerLineEdit.setText("NULL")


	def actualizePopMin(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.popMinLineEdit.text())
		try:
			if currentValue != "" and int(currentValue) >= 2 and int(currentValue) <= int(self.popMaxValue) :
				self.popMinValue = currentValue
			else:
				self.popMinValue = "1"
				self.ui.popMinLineEdit.setText("1")
		except Exception as e:
			self.displayError(error = str(e)+"\n Please enter an integer value !!!!")
			self.popMinValue = "1"
			self.ui.popMinLineEdit.setText("1")

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

	def actualizeRepMin(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.repMinLineEdit.text())
		try:
			if currentValue != "" and int(currentValue) >= 2 and int(currentValue) <= int(self.repMaxValue) :
				self.repMinValue = currentValue
			else:
				self.repMinValue = "1"
				self.ui.repMinLineEdit.setText("1")
		except Exception as e:
			self.displayError(error = str(e)+"\n Please enter an integer value !!!!")
			self.repMinValue = "1"
			self.ui.repMinLineEdit.setText("1")

	def actualizeRepMax(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.repMaxLineEdit.text())
		try:
			if currentValue != "" and int(currentValue) >= int(self.repMinValue):
				self.repMaxValue = currentValue
			else:
				self.repMaxValue = "10"
				self.ui.repMaxLineEdit.setText("10")
		except Exception as e:
			self.displayError(error = str(e)+"\n Please enter an integer value !!!!")
			self.repMaxValue = "10"
			self.ui.repMaxLineEdit.setText("10")

	def loadMatriceFile(self):
		"""Méthode qui permet de charger un fichier et afficher dans le plainText"""
		filename = QFileDialog.getOpenFileName(self, caption="Load the File", directory=os.getcwd(), filter="Text files (*.txt *.tab);;All (*.*)")
		if filename == "" and self.matricePathFile == "":
			self.displayError(error = "\"%s\" is not a valid file \n" % (filename))
		else:
			self.matricePathFile = filename
			self.workingDir = "/".join(str(filename).split("/")[:-1])+"/"
			self.basename = filename.split("/")[-1].split(".")[0]
			self.pathFileOut = "%s%s/" % (self.workingDir, self.basename)
			if os.path.isdir(self.pathFileOut):
				self.displayError(error = "Warnnig folder %s already exist,\nPlease remove/rename before run new analysis or use checkbox\n" % (self.pathFileOut))
			self.ui.matricePlainTextEdit.setPlainText(filename)
			self.actualizeRunBouton()

	def actualizeRunBouton(self):
		if self.matricePathFile != "" and self.nbIndiv != "NULL" and self.nbMarker != "NULL" :
			self.ui.runPushButton.setEnabled(True)
		else:
			self.ui.runPushButton.setDisabled(True)

	def resetPress(self):
		#"""reset PATH files"""
		self.ui.matricePlainTextEdit.clear()
		self.ui.frameRun.hide()
		self.ui.runningPlainTextEdit.clear()

		# réinisialise les champs
		self.ui.indivLineEdit.setText("NULL")
		self.ui.markerLineEdit.setText("NULL")
		self.ui.popMinLineEdit.setText("1")
		self.ui.popMaxLineEdit.setText("10")
		self.ui.repMinLineEdit.setText("1")
		self.ui.repMaxLineEdit.setText("10")

		self.ui.rmOldCheckBox.setChecked(False)

		# grise le bouton
		self.ui.runPushButton.setDisabled(True)

		# dégrise
		self.ui.loadMatriceFilePushButton.setEnabled(True)

		self.ui.indivLineEdit.setEnabled(True)
		self.ui.markerLineEdit.setEnabled(True)
		self.ui.popMinLineEdit.setEnabled(True)
		self.ui.popMaxLineEdit.setEnabled(True)
		self.ui.repMinLineEdit.setEnabled(True)
		self.ui.repMaxLineEdit.setEnabled(True)
		self.ui.rmOldCheckBox.setEnabled(True)
		self.ui.mainparamsPlainTextEdit.setEnabled(True)
		self.ui.extraparamsPlainTextEdit.setEnabled(True)

		# remet les variables a zero
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

			self.ui.indivLineEdit.setDisabled(True)
			self.ui.markerLineEdit.setDisabled(True)
			self.ui.popMinLineEdit.setDisabled(True)
			self.ui.popMaxLineEdit.setDisabled(True)
			self.ui.repMinLineEdit.setDisabled(True)
			self.ui.repMaxLineEdit.setDisabled(True)
			self.ui.rmOldCheckBox.setDisabled(True)
			self.ui.mainparamsPlainTextEdit.setDisabled(True)
			self.ui.extraparamsPlainTextEdit.setDisabled(True)


			self.mainparams = str(self.ui.mainparamsPlainTextEdit.toPlainText().toUtf8())
			self.extraparams = str(self.ui.extraparamsPlainTextEdit.toPlainText().toUtf8())

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

			# Construit les répertoires sh trash et .sge file
			outputSHDir = self.pathFileOut+"sh/"
			outputTrashDir = self.pathFileOut+"trash/"
			SGENameFile = self.pathFileOut+"submitQsubstructure.sge"

			txtInfo += """\
 // - Intput Info:
	- Input matrice is: %s"
	- Output prefix name is: %s
	- You want %s < K < %s and %s < Repetition < %s
 // - Output Info:
	- Working directory is: %s
	- Output directory is: %s\n
""" % (self.matricePathFile, self.basename, self.popMinValue, self.popMaxValue, self.repMinValue, self.repMaxValue, self.workingDir, self.pathFileOut)

			# création des répertoires et fichier mainparams
			os.makedirs(outputSHDir)											# création d'un dossier sh_scripts pour lancer les analyses structures
			os.makedirs(outputTrashDir)


			count=1
			for rep in range(int(self.repMinValue),int(self.repMaxValue)+1):																	# boucle sur le nombre de répétition
				os.makedirs(self.pathFileOut+"/repetition_"+str(rep))												# Création du répertoire correspondant
				for pop in range(int(self.popMinValue),int(self.popMaxValue)+1):																# boucle sur la variation du K (np pop)
					#print(str(pop))
					os.makedirs(self.pathFileOut+"/repetition_"+str(rep)+"/population_"+str(pop))					# création du répertoire de K
					mainparamsOut=open(self.pathFileOut+"/repetition_"+str(rep)+"/population_"+str(pop)+"/mainparams","w")		# ouverture du fichier mainparams correspondant
					extraparamsOut=open(self.pathFileOut+"/repetition_"+str(rep)+"/population_"+str(pop)+"/extraparams","w")		# ouverture du fichier extraparams correspondant (restera vide)

					# modifie mainarams pour adapter au numero de section , nsection et paths
					dictToReplace = {
					"**POP**"		:	str(pop),
					"**nIndiv**"	:	str(self.nbIndiv),
					"**nLoci**"		:	str(self.nbMarker),
					"**FILEIN**"	:	str(self.matricePathFile),
					"**FILEOUT**"	:	str(self.basename)+"_K"+str(pop)+"_run"+str(rep)+".txt"
					}

					mainparamsr = replace_all(dictToReplace, self.mainparams)

					mainparamsOut.write(mainparamsr)																			# Ecrit le mainparams

					mainparamsOut.close()

					extraparamsOut.write(self.extraparams)																			# Ecrit le mainparams
					extraparamsOut.close()
					#  écriture des scripts pour lancer les annalyses
					shOut=open(outputSHDir+"/"+str(count)+"_structure.sh","w")
					shOut.write("module load bioinfo/structure/2.3.4\n")
					shOut.write("cd "+self.pathFileOut+"/repetition_"+str(rep)+"/population_"+str(pop)+"/\n")
					shOut.write("structure")
					shOut.close()
					count+=1
					#  écriture du scripts qsub


			headerSGE = """\
#!/bin/bash

#$ -N structure
#$ -cwd
#$ -V
#$ -e """+outputTrashDir+"""
#$ -o """+outputTrashDir+"""
#$ -q long.q
#$ -t 1-"""+str(count-1)+"""
#$ -tc 100
#$ -S /bin/bash

/bin/bash """+outputSHDir+"""${SGE_TASK_ID}_structure.sh"""

			with open(SGENameFile, "w") as SGEFile:
				SGEFile.write(headerSGE)

			## Display a summary of the execution
			txtInfo +="""\
 // - Execution summary:
	- %s directories have been created corresponding to the number repeats,each with %s subdirectories corresponding to the variation number of populations (K).
	- To launch Structure execute:
qsub %s
	on the cluster.\n""" % (self.repMaxValue, self.popMaxValue, SGEFile.name)

			self.ui.statusbar.showMessage(str("FINISH, script product on : %s" % self.pathFileOut),9600)

			if args.cmdMode:
				txtInfo = txtInfo.replace("\t","\033[0m\t").replace("//","\033[92m").replace("qsub","\033[91m qsub")
				print(txtInfo)
			else:
				self.ui.runningPlainTextEdit.setPlainText(txtInfo)
				currentValue = str(self.ui.runningPlainTextEdit.toPlainText().toUtf8())
				highlight = syntax.PythonHighlighter(self.ui.runningPlainTextEdit.document())
				self.ui.runningPlainTextEdit.setPlainText(currentValue.decode("utf-8"))

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
	myapp = STRUCTURE(app)

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
	myapp = STRUCTURE(app)

	# Load info arguments to Class
	# matrice file
	myapp.matricePathFile = relativeToAbsolutePath(args.matriceParam)

	# info tab
	myapp.nbIndiv = args.nbIndivParam
	myapp.nbMarker = args.nbMarkerParam

	# K values
	myapp.popMinValue = args.nbpopiParam
	myapp.popMaxValue = args.nbpopmParam

	# rep values
	myapp.repMinValue = args.nbRepiParam
	myapp.repMaxValue = args.nbRepmParam

	# working dir path
	workingDir = "/".join(relativeToAbsolutePath(args.matriceParam).encode("utf-8").split("/")[:-1])+"/"
	myapp.workingDir = workingDir.encode("utf-8")
	# basename
	basename = relativeToAbsolutePath(args.matriceParam).encode("utf-8").split("/")[-1].split(".")[0]
	myapp.basename = basename.encode("utf-8")
	# pathFileOut dir path
	pathFileOut = workingDir+basename+"/"
	myapp.pathFileOut = pathFileOut.encode("utf-8")

	# rm old folder
	myapp.rmOld = args.rmOldParam

	# mainparams
	if args.mainFile == None:
		print("Use default mainparams\n")
		myapp.mainparams = mainparams
	else:
		with open(args.mainFile, "r") as fileO:
			main = "\n".join(fileO.readlines())
		myapp.mainparams = main

	# extraparams
	if args.extraFile == None:
		print("Use default extraparams\n")
		myapp.extraparams = extraparams
	else:
		with open(args.extraFile, "r") as fileO:
			extra = "\n".join(fileO.readlines())
		myapp.extraparams = extra

	# Run programme
	myapp.run()


if __name__ == '__main__':

	# Parameters recovery
	parser = argparse.ArgumentParser(prog='GUI_structure.py', description='''This programm build script to run multiple structure analysis.\n
																				#If use on cluster you can run in commande line with option -c and args''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display GUI_structure.py version number and exit')
	parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running if -c use')
	filesReq.add_argument('-c', '--cmd', action='store_true', dest = 'cmdMode', help = 'If used, programme run in CMD without interface')
	filesReq.add_argument('-i', '--matrice', metavar="<fileName>",type=existant_file, required=False, dest = 'matriceParam', help = 'input file matrice')
	filesReq.add_argument('-m', '--main', metavar="<fileName>",type=existant_file, required=False, dest = 'mainFile', help = 'mainparams file')
	filesReq.add_argument('-e', '--extra', metavar="<fileName>",type=existant_file, required=False, dest = 'extraFile', help = 'extraparams file')

	filesReq.add_argument('-nbi', '--nbIndiv', metavar="<int>",type = int, required=False, dest = 'nbIndivParam', help = 'Number of individus un matrice')
	filesReq.add_argument('-nbm', '--nbMarker', metavar="<int>",type = int, required=False, dest = 'nbMarkerParam', help = 'Number of markers un matrice')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-ri', '--repi', metavar="<int>",type = int, default=1, required=False, dest = 'nbRepiParam', help = 'Number of repetition min (default = 1)')
	files.add_argument('-rm', '--repm', metavar="<int>",type = int, default=10, required=False, dest = 'nbRepmParam', help = 'Number of repetition max (default = 10)')
	files.add_argument('-pi', '--popi', metavar="<int>",type = int, default=1, required=False, dest = 'nbpopiParam', help = 'Number of pop Min (default = 1)')
	files.add_argument('-pm', '--popm', metavar="<int>",type = int, default=10, required=False, dest = 'nbpopmParam', help = 'Number of pop Max  (default = 10)')
	files.add_argument('-o', '--outfile', metavar="<PrefixFileName>", required=False, dest = 'outputFile', help = 'output file Prefix (default = name of matrice file)')
	files.add_argument('-r', '--rm', metavar="<True/False>", choices=("True","False"), default="False", required=False, dest = 'rmOldParam', help = 'if directory exist remove (default = False)')	## Check parameters



	args = parser.parse_args()

	if args.cmdMode:
		for arg in [args.matriceParam, args.nbIndivParam, args.nbMarkerParam]:
			if arg == None :
				printCol.red("\nERROR: You must enter require arguments\n")
				parser.print_help()
				exit()
		cmd()
	else:
		# run interface
		main()
