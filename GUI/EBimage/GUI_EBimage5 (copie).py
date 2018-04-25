#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
## @package GUI_EBimage.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################



import argparse
import sys, os
import glob

# Python QT modules
from PyQt5.QtWidgets import QApplication, QFileSystemModel, QTableWidget, QFileDialog
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import QObject, pyqtSignal
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5 import uic

#import syntax
from rpy2.robjects.packages import importr
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

#exit()
##################################################
## Variables Globales
version = '1.0'
VERSION_DATE = '02/10/2017'

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
		self.ui.setWindowTitle("GUI_EBimage")
		self.ui.show()
		self.setFocusPolicy(QtCore.Qt.StrongFocus)
		self.activateWindow()
		self.initializationCalibrationVariables()
		self.initializationAnalysisVariables()
		self.initializationCommonVariables()
		self.createWidgets()


	def initializationCalibrationVariables(self):
		"""Initialize calibration variables with defaults values"""
		self.dicoFoldersCalibration =	{"leaf": None,
							 "symptom": None,
							 "background": None,
							}
							#"cat1": None,
							#"cat2": None,

		self.futherCat = 0
		self.dicoObjectOpenLineEditCalibration = {
						"leaf": self.ui.leafOpenLineEdit,
						"symptom": self.ui.symptomOpenLineEdit,
						"background": self.ui.backgroundOpenLineEdit,
					 }
						#"cat1": self.ui.category1OpenLineEdit,
						#"cat2": self.ui.category2OpenLineEdit
		self.dicoObjectOpenPushButtonCalibration = {
						"leaf": self.ui.loadLeafPushButton,
						"symptom": self.ui.loadSymptomPushButton,
						"background": self.ui.loadBackgroundPushButton,
					 }
						#"cat1": self.ui.loadCategory1PushButton,
						#"cat2": self.ui.loadCategoryPushButton
		self.calibrationOutPath = ""
		self.calibrationBasename = ""
		self.actualizeOutFiles()

	def initializationAnalysisVariables(self):
		"""Initialize analysis variables with defaults values"""
		self.leafMinSize = 1000
		self.leafBorderSize = 3
		self.symptomMinSize = 10
		self.symptomBorderSize = 3
		self.inputAnalyseFolder = ""
		self.outputAnalyseFolder = ""

	def initializationCommonVariables(self):
		"""Initialize common variables with defaults values"""
		#self.EBimageCalibration = EBimageCalibration
		#self.EBimageAnalysis = EBimageAnalysis
		self.tableWidget = QTableWidget()

	def actualizeOutFiles(self):
		"""actualise path of out files"""
		self.calibrationFilesOut = {
							"RData": self.calibrationOutPath+"/"+self.calibrationBasename+".RData",
							"png": self.calibrationOutPath+"/"+self.calibrationBasename+".png",
							"txt": self.calibrationOutPath+"/"+self.calibrationBasename+".txt",
							"info": self.calibrationOutPath+"/"+self.calibrationBasename+"_info.txt"
							}

	def createWidgets(self):
		"""Mise en place du masquage des frames et des connections de l'interface"""
		## Initialisation des frames et bouttons
		self.ui.runAnalysisPushButton.setDisabled(True)
		self.ui.runCalibrationPushButton.setDisabled(True)
		self.ui.runCalibrationFrame.hide()
		self.ui.calibrationOutFrame.hide()
		self.ui.futherFrame.hide()
		self.ui.cat1Frame.hide()
		self.ui.cat2Frame.hide()

		self.setWindowIcon(QIcon(resource_path("./includes/icon.ico")))
		self.ui.statusbar.setStyleSheet("color: rgb(255, 107, 8);font: 8pt 'Arial Black';")

		#self.resizeWindows()

		## Edition des connect:
		self.ui.loadCategory1PushButton.clicked.connect(lambda: self.loadFolder(inputCat = "cat1", methode="clicked"))
		self.ui.category1OpenLineEdit.editingFinished.connect(lambda: self.loadFolder(inputCat = "cat1", methode="write"))

		self.ui.loadCategoryPushButton.clicked.connect(lambda: self.loadFolder(inputCat = "cat2", methode="clicked"))
		self.ui.category2OpenLineEdit.editingFinished.connect(lambda: self.loadFolder(inputCat = "cat2", methode="write"))

		self.ui.loadLeafPushButton.clicked.connect(lambda: self.loadFolder(inputCat = "leaf", methode="clicked"))
		self.ui.leafOpenLineEdit.editingFinished.connect(lambda: self.loadFolder(inputCat = "leaf", methode="write"))

		self.ui.loadBackgroundPushButton.clicked.connect(lambda: self.loadFolder(inputCat = "background", methode="clicked"))
		self.ui.backgroundOpenLineEdit.editingFinished.connect(lambda: self.loadFolder(inputCat = "background", methode="write"))

		self.ui.loadSymptomPushButton.clicked.connect(lambda: self.loadFolder(inputCat = "symptom", methode="clicked"))
		self.ui.symptomOpenLineEdit.editingFinished.connect(lambda: self.loadFolder(inputCat = "symptom", methode="write"))

		self.ui.futherCategoriesComboBox.currentIndexChanged[int].connect(self.actualizeFutherCat)

		self.ui.runCalibrationPushButton.clicked.connect(self.run)
		self.ui.resetCalibrationPushButton.clicked.connect(self.resetLoadFolder)

	def resizeWindows(self):
		"""change la taille de fenetre"""
		# resize
		size = self.ui.centralwidget.sizeHint()
		self.ui.setGeometry(800,600,size.width(), size.height())

	def resetLoadFolder(self):
		"""To reset if delete of change value"""
		for inputCat in self.dicoFoldersCalibration.keys():
			self.dicoFoldersCalibration[inputCat] = None
			self.dicoObjectOpenLineEditCalibration[inputCat].setText("")
		self.actualizeRunButton()
		self.initializationCalibrationVariables()
		self.enableButtonsCalibration()
		self.ui.calibrationOutFrame.hide()
		self.viewLayout.removeWidget(self.tableWidget)

	def loadFolder(self,inputCat = None, methode = None):
		"""Méthode qui permet de charger un fichier et afficher dans le plainText"""
		directoryToOpen = os.getcwd()
		for key, folder in self.dicoFoldersCalibration.items():
			if folder != None:
				#print(key, folder)
				directoryToOpen = folder
				break

		if methode == "clicked":
			pathdir = QFileDialog.getExistingDirectory(self, caption="Load the directory "+inputCat, directory=directoryToOpen)
			self.dicoObjectOpenLineEditCalibration[inputCat].setText(pathdir)
		elif methode == "write":
			pathdir = str(self.dicoObjectOpenLineEditCalibration[inputCat].text())

		if pathdir != "" and os.path.isdir(pathdir):
			self.dicoFoldersCalibration[inputCat] = pathdir

		else:
			self.dicoFoldersCalibration[inputCat] = ""
			self.dicoObjectOpenLineEditCalibration[inputCat].setText("")
			self.displayError(error = "\"%s\" is not a valid Path \n" % (pathdir))
		self.actualizeRunButton()
		self.autocomplete_Folder()
		#print(self.dicoFoldersCalibration[inputCat])

	def autocomplete_Folder(self):
		"""autocomplet les répertoires si un déja remplit"""
		directoryToOpen = None
		for key, folder in self.dicoFoldersCalibration.items():
			if folder != None:
				directoryToOpen = "/".join(str(folder).split("/")[:-1])
				break
		if directoryToOpen != None:
			for root, dirs, files in os.walk(directoryToOpen):
				#if args.debug: print(root, dirs, files)
				for inputCat in self.dicoFoldersCalibration.keys():
					if inputCat in dirs:
						pathdir = "{}/{}".format(root,inputCat)
						self.dicoFoldersCalibration[inputCat] = pathdir
						self.dicoObjectOpenLineEditCalibration[inputCat].setText(pathdir)

		print(directoryToOpen)
		# répertoire racine pour le QTreeView
		rephome = directoryToOpen

		#-- Modèle
		self.myModel = QFileSystemModel()
		self.myModel.setReadOnly(False)
		self.myModel.setRootPath(rephome)

		#-- treeview
		self.treeView.setModel(self.myModel)

		rootModelIndex = self.myModel.index(rephome)
		self.treeView.setRootIndex(rootModelIndex)

















		self.actualizeRunButton()


	def actualizeFutherCat(self):
		"""masque les autres catégorie quand changer"""
		oldState = self.futherCat
		self.futherCat  = str(self.ui.futherCategoriesComboBox.currentText())
		if self.futherCat == "0":
			self.ui.cat1Frame.hide()
			self.ui.cat2Frame.hide()
			self.dicoFoldersCalibration["cat1"] = None
			self.dicoFoldersCalibration["cat2"] = None
			self.dicoObjectOpenLineEditCalibration["cat1"].setText("")
			self.dicoObjectOpenLineEditCalibration["cat2"].setText("")
			self.displayError(error = "Reset Futher categories path !! \n")
		elif self.futherCat == "1":
			self.ui.cat1Frame.show()
			self.ui.cat2Frame.hide()
			if oldState == "2":
				self.dicoObjectOpenLineEditCalibration["cat2"].setText("")
				self.displayError(error = "Reset categories 2 path !! \n")
			self.dicoFoldersCalibration["cat2"] = None
		elif self.futherCat == "2":
			self.ui.cat1Frame.show()
			self.ui.cat2Frame.show()
		self.actualizeRunButton()



	def actualizeRunButton(self):
		"""de grise le bouton run si rempli matrice et order"""
		if self.futherCat == 0:
			for key, folder in self.dicoFoldersCalibration.items():
				#print(key, folder)
				if folder != None :
					self.ui.runCalibrationPushButton.setEnabled(True)
				elif key not in ["cat1", "cat2"]:
					self.ui.runCalibrationPushButton.setDisabled(True)
					break
		else:
			if self.futherCat == "1" and self.dicoFoldersCalibration["cat1"] != None:
				self.ui.runCalibrationPushButton.setEnabled(True)
			elif self.futherCat == "2" and self.dicoFoldersCalibration["cat1"] != None and self.dicoFoldersCalibration["cat2"] != None:
				self.ui.runCalibrationPushButton.setEnabled(True)
			else:
				self.ui.runCalibrationPushButton.setDisabled(True)


	def infoDialogue(self, status): ## Method to open a message box
			infoBox = QMessageBox() ##Message Box that doesn't run
			infoBox.setFixedSize(1000,1000)
			infoBox.setIcon(QMessageBox.Information)
			infoBox.setText("Calibration OK")
			if status == "new":
				infoBox.setInformativeText("all files created")
			if status == "already":
				infoBox.setInformativeText("all files already exist")
			infoBox.setWindowTitle("Good Calibration")
			filestxt = "Files created:\n"
			for key, path in self.calibrationFilesOut.items():
				filestxt += "- {}\n".format(path)
			infoBox.setDetailedText(filestxt)
			infoBox.setStandardButtons(QMessageBox.Ok)
			infoBox.setEscapeButton(QMessageBox.Close)
			infoBox.exec_()
			self.openCalibrationTable()
			self.loadCalibrationPNG()

	def disalbledButtonsCalibration(self):
		"""Disable buttons when run calibration"""
		self.ui.runCalibrationPushButton.setDisabled(True)
		for values in self.dicoObjectOpenPushButtonCalibration.values():
			values.setDisabled(True)
		for values in self.dicoObjectOpenLineEditCalibration.values():
			values.setDisabled(True)

	def enableButtonsCalibration(self):
		"""enable buttons when finish run calibration or already do"""
		self.ui.runCalibrationPushButton.setDisabled(True)
		for values in self.dicoObjectOpenPushButtonCalibration.values():
			values.setEnabled(True)
		for values in self.dicoObjectOpenLineEditCalibration.values():
			values.setEnabled(True)
			self.ui.resetCalibrationPushButton.setEnabled(True)

	def run(self):
		try:
			warning = ""	# initialise le nombre d'erreur
			val = 0			# initialise le nombre d'erreur
			txtInfo = ""
			result = "0"
			# grise
			self.disalbledButtonsCalibration()
			self.ui.resetCalibrationPushButton.setDisabled(True)

			# import R's "base" package
			base = importr('base')

			with open("fonctions_apprentissage.r", "r", encoding="utf-8") as apprentissageRopen:
				apprentissage = "".join(apprentissageRopen.readlines())

			apprentissage = SignatureTranslatedAnonymousPackage(apprentissage, "apprentissage")

			self.calibrationOutPath = "/".join(str(self.dicoFoldersCalibration["leaf"]).split("/")[:-1])
			self.calibrationBasename = self.calibrationOutPath.split("/")[-1]
			self.actualizeOutFiles()

			if args.debug: print("{}\n{}".format(apprentissage, dir(apprentissage)))

			# test if Rdata file already exist, if yes remove file if user say yes, or stop analyse
			if os.path.exists(self.calibrationFilesOut["RData"]):
				reply = QMessageBox.question(self, 'Warning', 'File will be overwritten.\nDo you still want to proceed?', QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
				if reply == QMessageBox.Yes:
					for key, path in self.calibrationFilesOut.items():
						os.remove(path)
					reloadCalibration = True
				elif reply == QMessageBox.No:
					reloadCalibration = False
			else:
				reloadCalibration = True

			if reloadCalibration:
				self.ui.statusbar.showMessage(str("Running calibration, please waiting ...."),9600)
				#result , self.calibrationFilesOut["RData"] = apprentissage.apprentissage(self.dicoObjectOpenLineEditCalibration["leaf"],self.dicoObjectOpenLineEditCalibration["symptom"],self.dicoObjectOpenLineEditCalibration["background"]).r_repr().replace('"','').replace("c(","").replace(")","").split(",")
				#doivent être (dans cet ordre) le nom (relatif) des sous-répertoires fond, limbe, lésions
				backName, leafName, symptomName = self.dicoFoldersCalibration["background"].split("/")[-1], self.dicoFoldersCalibration["leaf"].split("/")[-1], self.dicoFoldersCalibration["symptom"].split("/")[-1]
				result , good = apprentissage.apprentissage(self.calibrationOutPath,backName,leafName,symptomName).r_repr().replace('"','').replace("c(","").replace(")","").split(",")
				self.calibrationFileOpenLineEdit.setText(self.calibrationFilesOut["RData"])
			if result == "1" and os.path.exists(self.calibrationFilesOut["RData"]):
				print(result, self.calibrationFilesOut["RData"])
				self.infoDialogue(status = "new")
				self.ui.statusbar.showMessage(str("FINISH, files were product on : %s" % self.calibrationOutPath),9600)
				self.ui.resetCalibrationPushButton.setEnabled(True)
			elif  result == "0" and os.path.exists(self.calibrationFilesOut["RData"]):
				self.infoDialogue(status = "already")
				print(result, self.calibrationFilesOut["RData"])
				self.calibrationFileOpenLineEdit.setText("")
				self.ui.resetCalibrationPushButton.setEnabled(True)
				self.resetLoadFolder()
				self.enableButtonsCalibration()
			elif result == "0" and not os.path.exists(self.calibrationFilesOut["RData"]):
				self.displayError(error = "Error when running R code....")

		except Exception as e:
			self.displayError(error = e)

	def openCalibrationTable(self):
		"""View result of calibration"""
		self.ui.calibrationOutFrame.show()
		with open(self.calibrationFilesOut["info"]) as tableFile:
			lines = [line.rstrip().split() for line in tableFile.readlines()]
		self.tableWidget.setRowCount(len(lines)-1)
		self.tableWidget.setColumnCount(4)
		row = 0
		print(lines[1:])
		for line in lines[1:]:
			print(line)
			if len(line) == 0:
				next
			if len(line) == 3:
				self.tableWidget.setItem(row,0, QTableWidgetItem(""))
				self.tableWidget.setItem(row,1, QTableWidgetItem(lines[row+1][0]))
				self.tableWidget.setItem(row,2, QTableWidgetItem(lines[row+1][1]))
				self.tableWidget.setItem(row,3, QTableWidgetItem(lines[row+1][2]))
			if len(line) == 4:
				self.tableWidget.setItem(row,0, QTableWidgetItem(lines[row+1][0]))
				self.tableWidget.setItem(row,1, QTableWidgetItem(lines[row+1][1]))
				self.tableWidget.setItem(row,2, QTableWidgetItem(lines[row+1][2]))
				self.tableWidget.setItem(row,3, QTableWidgetItem(lines[row+1][3]))
			row += 1
		self.viewLayout.addWidget(self.tableWidget)
		self.show()
		self.loadCalibrationPNG()

	def loadCalibrationPNG(self):
		"""Load png file if good calibration"""
		print(self.calibrationFilesOut["png"])
		pic = QPixmap(self.calibrationFilesOut["png"])
		pic = pic.scaled(480,480)
		self.ui.imgLabel.setPixmap(pic)


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
	parser.add_argument('-d', '--debug',action='store_true', help='enter verbose/debug mode', dest = "debug", default = "False")

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
