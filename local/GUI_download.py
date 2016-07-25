#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
## @package buildTableauResume.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################

import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath, directory, existant_file, loadInList

# Python modules
import argparse
import os, sys,shutil
import re  # Regular expression library
import urllib

from subprocess import Popen

from Bio import SeqIO
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import  *
from PyQt4.QtGui import *
from PyQt4 import uic


##################################################
## Variables Globales
version='1.0'
VERSION_DATE='25/07/2016'
_toolname = 'genbank-download'
_email = 'dev@simon.net.nz'

def resource_path(relative):
	if hasattr(sys, "_MEIPASS"):
		return os.path.join(sys._MEIPASS, relative)
	return os.path.join(relative)



############################ FOR MAC PARM ######################################
# les .app démarrent python avec '/' comme cwd
if "darwin" in sys.platform : # and ".app/" in sys.argv[0]:
	# pour aller a l'interieur du .app
	mycwd = sys.argv[0].split(".app/")[0] + ".app/Contents/Resources/"
	os.chdir(mycwd)


class IDError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)


filename = 'download.ui'
myfontfile = resource_path(os.path.join(filename))
formDownload,baseDownload = uic.loadUiType(myfontfile)



# ************************************* CLASSE Download Gestion de l'affichage et de la récupération de donnée ****************************************************
## @class Download
# @brief Classe principale, fenêtre principale
class Download( formDownload, baseDownload ):
	""" Classe principale qui est aussi la fenêtre principale de la GUI
	"""
	def __init__(self,app,parent=None):
		super(Download,self).__init__(parent)
		self.app = app
		self.initialiseVariables()
		self.createWidgets()

	def initialiseVariables(self):
		"""Initialise les variables par défaults"""
		self.pathFileOut = "./"
		self.pathFileOutGB = "Genbank_Files"
		self.outputFileNucl = "Nucleique.fasta"
		self.outputFileProt = "Proteique.fasta"
		self.listID = []
		self.listIDextract = []
		self.typeChoice = "Both Nucl and Prot"
		self.keepCheck = "False"
		self.buildCheck = "False"
		# pour les fichiers de sortie
		self.listRecordsProt = []
		self.listRecordsNucl = []
		self.tableau = []


	def createWidgets(self):
		"""Mise en place du masquage des frames et des connections de l'interface"""
		self.ui = self
		self.ui.setupUi(self)
		## Initialisation des frames et bouttons
		self.ui.displayErrorEdit.hide()
		self.ui.progressbar.hide()


		## Edition des connect:
		QObject.connect(self.ui.loadFilePushButton,SIGNAL("clicked()"),self.loadFile)
		QObject.connect(self.ui.resetPushButton,SIGNAL("clicked()"),self.resetPress)
		QObject.connect(self.ui.runPushButton,SIGNAL("clicked()"),self.run)

		QObject.connect(self.ui.choiceComboBox,SIGNAL("currentIndexChanged(int)"),self.actualizeChoice)

		QObject.connect(self.ui.keepGBcheckBox,SIGNAL("stateChanged(int)"),self.actualizeKeep)
		QObject.connect(self.ui.buildTabcheckBox,SIGNAL("stateChanged(int)"),self.actualizeBuildTab)

		QObject.connect(self.ui.nuclOutLineEdit,SIGNAL("editingFinished()"),self.actualizeOutputNucl)
		QObject.connect(self.ui.protOutLineEdit,SIGNAL("editingFinished()"),self.actualizeOutputProt)
		QObject.connect(self.ui.pathOutlineEdit,SIGNAL("editingFinished()"),self.actualizepathOutline)


	def actualizeKeep(self):
		"""change la valeur du choix quand changer"""
		if self.ui.keepGBcheckBox.isChecked():
			self.keepCheck = "True"
		else:
			self.keepCheck = "False"
	def actualizeBuildTab(self):
		"""change la valeur du choix quand changer"""
		if self.ui.buildTabcheckBox.isChecked():
			self.buildCheck = "True"
		else:
			self.buildCheck = "False"

	def actualizeChoice(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.choiceComboBox.currentText())
		self.typeChoice = currentValue
		if self.typeChoice in ["Nucleique Only"]:
			self.ui.nuclOutLineEdit.setEnabled(True)
			self.ui.protOutLineEdit.setDisabled(True)
		elif self.typeChoice in ["Proteique Only"]:
			self.ui.nuclOutLineEdit.setDisabled(True)
			self.ui.protOutLineEdit.setEnabled(True)
		else:
			self.ui.nuclOutLineEdit.setEnabled(True)
			self.ui.protOutLineEdit.setEnabled(True)

	def actualizeOutputNucl(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.nuclOutLineEdit.text())
		if currentValue != "":
			self.outputFileNucl = currentValue
		else:
			self.outputFileNucl = "Nucleique.fasta"
			self.ui.nuclOutLineEdit.setText("Nucleique.fasta")

	def actualizeOutputProt(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.protOutLineEdit.text())
		if currentValue != "":
			self.outputFileProt = currentValue
		else:
			self.outputFileProt = "Proteique.fasta"
			self.ui.protOutLineEdit.setText("Proteique.fasta")

	def actualizepathOutline(self):
		"""change la valeur du choix quand changer"""
		currentValue = str(self.ui.pathOutlineEdit.text())
		if currentValue != "":
			self.pathFileOut = currentValue
		else:
			self.pathFileOut = "./"
			self.ui.pathOutlineEdit.setText("./")


	def string_cleanup(self, x, notwanted, replaceFor):
		"""Sert a cleaner les txt"""
		for item in notwanted:
			x = re.sub(item, replaceFor, x)
		return x

	def formatTxt(self, txtload):
		"""Sert a formater les donnees pour mettre les id en liste"""
		# replace les caracteres espace et tab par des retours a la ligne
		txtformat = self.string_cleanup(txtload, ["\s+", "\t", "\r+"], "\n")
		txtformat2 = self.string_cleanup(txtformat, ["\n{2,}", "\r{2,}"], "\n").rstrip("\r\n")
		return txtformat2


	def loadFile(self):
		"""Méthode qui permet de charger un fichier et afficher dans le plainText"""
		filename = QtGui.QFileDialog.getOpenFileName(self, caption="Load the File", directory=os.getcwd(), filter="Text file *.txt")
		if filename == "":

			self.displayErrorLoad(error = "\"%s\" is not a valid file \n" % (filename))
		else:
			list = open(filename,"r").readlines()
			listgood=[line.rstrip() for line in list]
			self.ui.statusbar.showMessage(str("The file is load from %s" % filename),7200)
			txtload = "\n".join(listgood)
			txtformat = self.formatTxt(txtload)
			self.ui.IDplainTextEdit.setPlainText(txtformat)


	def resetPress(self):
		"""reset load or list ID"""
		self.ui.IDplainTextEdit.clear()
		self.ui.nuclOutLineEdit.setText("Nucleique.fasta")
		self.ui.protOutLineEdit.setText("Proteique.fasta")
		self.ui.pathOutlineEdit.setText("./")
		self.ui.progressbar.hide()
		self.ui.runPushButton.setEnabled(True)
		self.ui.displayErrorEdit.hide()
		self.ui.displayErrorEdit.clear()

		# uncheck les checkBoxs
		self.ui.keepGBcheckBox.setChecked(False)
		self.ui.buildTabcheckBox.setChecked(False)

		self.ui.runPushButton.setEnabled(True)
		self.ui.IDplainTextEdit.setEnabled(True)
		self.ui.nuclOutLineEdit.setEnabled(True)
		self.ui.protOutLineEdit.setEnabled(True)
		self.ui.pathOutlineEdit.setEnabled(True)
		self.ui.choiceComboBox.setEnabled(True)
		self.ui.keepGBcheckBox.setEnabled(True)
		self.ui.buildTabcheckBox.setEnabled(True)
		self.ui.loadFilePushButton.setEnabled(True)
		#if os.path.isdir(self.pathFileOutGB):
			#shutil.rmtree(self.pathFileOutGB)
		self.initialiseVariables()



	def get_accession(self,query, database, rettype):
		"""
		Returns a nucleotide sequence from genbank.

		:param query: the accession number
		:param database: the database to use (default=nucleotide)
		:param rettype: the return type for the sequence (native,fasta,gb,xml)

		:return: text of sequence in requested `rettype`
		:rtype: string

		"""

		params = {
			'db': database,
			'tool': _toolname,
			'email': _email,
			'id': query,
			'rettype': rettype,
		}
		url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
		url = url + urllib.urlencode(params)
		try:
			#print url
			data = urllib.urlopen(url).read()
			return data
		except Exception as erreurs:
			print erreurs
			raise IDError("\"%s\" is not a valid accession (error:%s) \n" % (query, erreurs))

	def removeGB(self):
		"""For remove files GB if not check keep"""
		if self.keepCheck == "False":
			shutil.rmtree(self.pathFileOutGB)

	def txtToGenbankObject(self,gbId,geneBankFile):
		tempGBfile = open(self.pathFileOutGB+gbId, "w")
		tempGBfile.write(geneBankFile)
		tempGBfile.close()
		tempGBfile2 = open(self.pathFileOutGB+gbId, "r")
		gb=SeqIO.parse(tempGBfile2,'genbank') #you MUST tell SeqIO what format is being read)

		return gb

	def run(self):
		try:
			warning = ""	# initialise le nombre d'erreur
			val = 0			# initialise le nombre d'erreur
			if args.cmdMode:
				#print self.listID
				os.mkdir(self.pathFileOutGB)

			else:
				"""to run programme"""
				currentValue = str(self.ui.IDplainTextEdit.toPlainText())
				txtformat = self.formatTxt(currentValue)
				self.listID = txtformat.split("\n")


				# création du pathout GB
				if self.pathFileOut[-1] not in ["/","\\"]:
					if "darwin" in sys.platform :											# Pour l'utilisation sous Mac OS:
						self.pathFileOut += "/"
						self.pathFileOutGB = self.pathFileOut+"/Genbank_Files/"
					if "linux" in sys.platform:												# Pour l'utilisation sous Linux:
						self.pathFileOut += "/"
						self.pathFileOutGB = self.pathFileOut+"/Genbank_Files/"
					if "win" in sys.platform and "darwin" not in sys.platform:				# Pour l'utilisation sous windows:
						self.pathFileOut += "\\"
						self.pathFileOutGB = self.pathFileOut+"\\Genbank_Files\\"
				else:
					if "darwin" in sys.platform :											# Pour l'utilisation sous Mac OS:
						self.pathFileOutGB = self.pathFileOut+"/Genbank_Files/"
					if "linux" in sys.platform:												# Pour l'utilisation sous Linux:
						self.pathFileOutGB = self.pathFileOut+"/Genbank_Files/"
					if "win" in sys.platform and "darwin" not in sys.platform:				# Pour l'utilisation sous windows:
						self.pathFileOutGB = self.pathFileOut+"\\Genbank_Files\\"
				if os.path.isdir(self.pathFileOutGB):
					#print("Warnnig folder "+self.pathFileOutGB+" already exist")
					warning += "Warnnig folder "+self.pathFileOutGB+" already exist "
				else:
					os.mkdir(self.pathFileOutGB)

				#grise le bouton pour pas relancer job
				self.ui.runPushButton.setDisabled(True)
				self.ui.IDplainTextEdit.setDisabled(True)
				self.ui.nuclOutLineEdit.setDisabled(True)
				self.ui.protOutLineEdit.setDisabled(True)
				self.ui.pathOutlineEdit.setDisabled(True)
				self.ui.choiceComboBox.setDisabled(True)
				self.ui.keepGBcheckBox.setDisabled(True)
				self.ui.buildTabcheckBox.setDisabled(True)
				self.ui.loadFilePushButton.setDisabled(True)
				self.ui.displayErrorEdit.hide()
				self.ui.displayErrorEdit.clear()
				self.ui.progressbar.show()
				self.progressbar.setRange(0,len(self.listID)-1)
				self.progressbar.setValue(val)

				#self.ui.statusbar.showMessage(str("running programm with %i ID" % len(self.listID)-1),7200)

			# parcours de la liste des ID à chercher:
			for gbId in self.listID:
				genbankFicheID = self.get_accession(gbId, "nucleotide", "gb").decode("utf-8")		# Récupère la fiche genbank assicier à l'ID en txt
				if "empty" in genbankFicheID:
					#self.displayError(error = "\"%s\" is not a valid accession \n" % (gbId))
					raise IDError("\"%s\" is not a valid accession \n" % (gbId))

				gbfiche = self.txtToGenbankObject(gbId,genbankFicheID)							# convertie le txt en objet genebank

				# parcours la fiche genbank
				for record in gbfiche :
					val+=1
					seq = record.seq
					type = record.seq.alphabet

					if self.ui.buildTabcheckBox.isChecked() or args.paramtaxafile:
						################## POUR Denis #################
						taxonomy = "; ".join(record.annotations["taxonomy"])
						description = record.description
						organisme = record.annotations["organism"]
						accession = record.id
						txtTab = accession+"\t"+organisme+"\t"+taxonomy+"\t"+description
						#print(txtTab)
						self.tableau.append(txtTab)
						################## POUR Denis #################
					else:

						if str(type) == "IUPACProtein()": 											# La sequence est proteique on récupère l'accession du CDS
							self.listRecordsProt.append(record)										# Sauvegarde la fiche dans la liste prot
							features=record.features
							t=len(features)-1
							if "coded_by"  in features[t].qualifiers.keys():
								coded_by=features[t].qualifiers["coded_by"][0].split(":")[0].replace("complement(","").replace("join(","")
								if coded_by == "":
									coded_by=feature.qualifiers["coded_by"][0]
								genbankFicheIDNucl = self.get_accession(coded_by, "nucleotide", "gb").decode("utf-8")		# Récupère la fiche genbank assicier à l'ID en txt
								gbnucle = self.txtToGenbankObject(coded_by,genbankFicheIDNucl)							# convertie le txt en objet genebank
								for record in gbnucle:
									self.listRecordsNucl.append(record)										# Sauvegarde la fiche dans la liste Nucle
							#else:
								#print("WARNING: No coded_by for ID: %s" % gbId)
								#warning += "WARNING: No coded_by for ID: %s\n" % gbId

						else:																	# La séquence est nucleique on récupère l'accession prot
							self.listRecordsNucl.append(record)										# Sauvegarde la fiche dans la liste Nucle
							features=record.features
							t=len(features)-1
							for feature in features:
								if "protein_id"  in feature.qualifiers.keys():
									protein_id=feature.qualifiers["protein_id"][0].split(":")[0].replace("complement(","").replace("join(","")
									if protein_id == "":
										protein_id=feature.qualifiers["protein_id"][0]
									genbankFicheIDprot = self.get_accession(protein_id, "protein", "gb").decode("utf-8")		# Récupère la fiche genbank assicier à l'ID en txt
									gbprot = self.txtToGenbankObject(protein_id,genbankFicheIDprot)							# convertie le txt en objet genebank
									for record in gbprot:
										self.listRecordsProt.append(record)												# Sauvegarde la fiche dans la liste prot
								#else:
									#print("WARNING: No protein_id for ID: %s" % gbId)
									#warning += "WARNING: No protein_id for ID: %s\n" % gbId
					self.progressbar.setValue(val)
					self.listIDextract.append(gbId)
				gbfiche.close()
							# si des erreurs:
				if len(warning)>0:
					self.displayErrorLoad(error = warning)
			if self.ui.buildTabcheckBox.isChecked() or args.paramtaxafile:
				self.buildTabfile()
			else:
				self.printToFiles()

		except IDError as e:
			print e
			self.displayError(error = e)

	def printToFiles(self):
		"""Ajoute les téléchargement dans les fichiers"""

		if self.pathFileOut[-1] not in ["/","\\"]:
			# Pour l'utilisation sous Mac OS:
			if "darwin" in sys.platform :
				self.pathFileOut += "/"
			# Pour l'utilisation sous Linux:
			if "linux" in sys.platform:
				self.pathFileOut += "/"
			# Pour l'utilisation sous windows:
			if "win" in sys.platform and "darwin" not in sys.platform:
				self.pathFileOut += "\\"

		if self.typeChoice in ["Nucleique Only","Both Nucl and Prot"]:
			nbSeqNucle=0
			output_handleNucl = open(self.pathFileOut+self.outputFileNucl, "w")
			for recordNucl in self.listRecordsNucl:
				nbSeqNucle+=1
				SeqIO.write(recordNucl, output_handleNucl, "fasta")
			output_handleNucl.close()

		if self.typeChoice in ["Proteique Only","Both Nucl and Prot"]:
			nbSeqProt=0
			output_handleProt = open(self.pathFileOut+self.outputFileProt, "w")
			for recordProt in self.listRecordsProt:
				nbSeqProt+=1
				SeqIO.write(recordProt, output_handleProt, "fasta")
			output_handleProt.close()
		self.removeGB()

	def buildTabfile(self):
		"""build tab file if checkbox"""
		if self.buildCheck == "True":
			output_tabFile = open(self.pathFileOut+"tabTaxonomybuild.tab", "w")
			txtOut = "\n".join(self.tableau)
			output_tabFile.write(txtOut)
			output_tabFile.close()
		self.removeGB()


	def displayErrorLoad(self, error):
		""" affiche les erreurs dans la zone de text"""
		self.ui.displayErrorEdit.show()
		self.ui.displayErrorEdit.setPlainText(error)

	def displayError(self, error):
		""" affiche les erreurs dans la zone de text"""
		self.printToFiles()
		# Grise les cases pour ne pas relancer dessus et faire un reset
		self.ui.runPushButton.setDisabled(True)
		self.ui.IDplainTextEdit.setDisabled(True)
		self.ui.nuclOutLineEdit.setDisabled(True)
		self.ui.protOutLineEdit.setDisabled(True)
		self.ui.pathOutlineEdit.setDisabled(True)
		self.ui.choiceComboBox.setDisabled(True)
		self.ui.loadFilePushButton.setDisabled(True)
		self.ui.keepGBcheckBox.setDisabled(True)
		self.ui.buildTabcheckBox.setDisabled(True)

		self.ui.displayErrorEdit.show()
		txtIDPass = "Only this accession was download:\n%s" % "\n".join(self.listIDextract)
		txtError = error.value +"\n"+txtIDPass
		self.ui.displayErrorEdit.setPlainText(txtError)





def main():

	#print sys.argv+["-cmd", "-gnome-terminal"]
	nargv = sys.argv

	# instanciation des objets principaux
	app = QApplication(nargv)
	myapp = Download(app)

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
	myapp = Download(app)

	# Load info arguments to Class
	# List ID
	listIDLoad = loadInList(args.listFile)
	myapp.listID = listIDLoad
	#OUT PATH
	pathOut = relativeToAbsolutePath(args.paramOutPath)+"/"
	myapp.pathFileOut = pathOut
	#GB PATH
	myapp.pathFileOutGB = pathOut+"Genbank_Files/"
	# keep GB files:
	if args.paramkeepfile:
		myapp.keepCheck = "True"
	# build taxonomy file
	if args.paramtaxafile:
		myapp.buildCheck = "True"
	# build files choices
	if args.parambuildfile == "B":
		myapp.typeChoice = "Both Nucl and Prot"
	if args.parambuildfile == "N":
		myapp.typeChoice = "Nucleique Only"
	if args.parambuildfile == "P":
		myapp.typeChoice = "Proteique Only"

	# Run programme
	myapp.run()


if __name__ == '__main__':

	# Parameters recovery
	parser = argparse.ArgumentParser(prog='concatFastasFile.py', description='''This Programme open GUI to download genbank file with list of accession.\n
																				You can both build a taxonomy tab with ID.\n
																				If use on cluster you can run in cammande line with option -c and args''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display concatFastasFile.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running if -c use')
	filesReq.add_argument('-c', '--cmd', action='store_true', dest = 'cmdMode', help = 'If used, programme run in CMD without interface')
	filesReq.add_argument('-l', '--list', metavar="<filename>",type=existant_file, required=False, dest = 'listFile', help = 'List of ID')
	filesReq.add_argument('-o', '--out', metavar="<path/to/directory>",type=directory, required=False, dest = 'paramOutPath', help = 'Name of output directory')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-k', '--keepGB', action='store_true', required=False, dest = 'paramkeepfile', help = 'keepGBfile if add')
	files.add_argument('-t', '--taxa', action='store_true', required=False, dest = 'paramtaxafile', help = 'build taxonomy file if add')
	files.add_argument('-b', '--build', metavar="<B;N;P>", choices=("B","N", "P"),default="B", required=False, dest = 'parambuildfile', help = 'build Both (B) default, Nucleique (N) or Proteique (P) files')
	# Check parameters
	args = parser.parse_args()

	if args.cmdMode:
		cmd()
	else:
		# run interface
		main()

