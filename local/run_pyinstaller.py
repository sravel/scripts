#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
# @package run_pyinstaller.py
# @author Sebastien Ravel


"""
	The run_pyinstaller script
	==========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 15/09/2016
	:version: 0.1

	Script description
	------------------

	This run pyinstaller to build file on OS

	Example
	-------

	>>> run_pyinstaller.py -p GUI_DAPC.py -u ./includes/dapc.ui -a includes/ -i includes/icon.ico

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display run_raxml.py version number and exit

	Input mandatory infos for running:
		- \-p <filename>, --python <filename>
						Main programm script
		- \-u <filename>, --ui <filename>
						UI file associated to python
		- \-a <path/to/directory>, --addIncludes <path/to/directory>
						Name of directory with includes file (must be name
						includes)
		- \-i <filename>, --icon <filename>
						icon.ico file

	Input infos for running with default values:
		- \-o <path/to/directory>, --out <path/to/directory>
						Name of output file directory

"""


##################################################
## Modules
##################################################

## Python modules
import argparse, glob
from time import localtime, strftime
import sys, os, shutil, re

class directory(str):
	"""
	directory() CLASS
	=================

	Class which derives from string.
	Checks that the string is and path to valid directory and not empty

	Returns object which able to return basename and list file (with an extention, else consider as folder)

	Example:

	>>> inDirectory=directory("/home/sravel/Documents")
	>>> print(inDirectory.pathDirectory())
	>>> /home/sravel/Documents/

	>>> print(inDirectory.listFiles())
	>>> ["File1.txt","File2.pl","file.toto"]
	"""
	def __init__(self, pathDirectory = None):
		"""
			Initialise variable
		"""
		self.listPath = []			# all in the path
		self.listDir = []			# only directory in path
		self.listFiles = []			# only files in path

		self.current_dir = os.path.dirname(os.path.abspath(__file__))

		#Change relative path to absolute path
		self.pathDirectory = relativeToAbsolutePath(pathDirectory)

		# appel les fonctions
		self.testDirExist()
		self.lsInDir()
		self.splitFilesDir()

	def __repr__(self):
		"""Fonction qui permet de formater le text de sortie lors du print du dictionnaire"""
		txtOut = """
pathDirectory=%s\n
listPath=%s\n
listDir=%s\n
listFiles=%s\n
""" % (self.pathDirectory, str(self.listPath), str(self.listDir), str(self.listFiles))
		return txtOut

	def testDirExist(self):
		"""Test l'existance du répertoire"""
		if os.path.isdir(self.pathDirectory) != True :
			print("ERROR MODULES_SEB::Class-directory : path '%s' is not valide path" % self.pathDirectory )
			exit()

	def lsInDir(self):
		"""List all in directory"""
		self.testDirExist()
		if self.pathDirectory[-1] != "/":
			self.pathDirectory += "/"
		if self.pathDirectory[-1] != "*":
			pathDirectoryList = self.pathDirectory+"*"
		self.listPath=glob.glob(pathDirectoryList)

	def lsExtInDirToList(self, ext):
		"""List specific extention file in directory"""
		lsFilesFasta =[]
		for fichier in self.listPath:
			try:
				if "." in fichier.split("/")[-1]:
					nameFichier = fichier.split("/")[-1].split(".")[0]
					extentionFichier = fichier.split("/")[-1].split(".")[-1]
				else:
					nameFichier = fichier.split("/")[-1]
					extentionFichier = ""
			except:
				extentionFichier = "directory"
			if extentionFichier in ext:
				lsFilesFasta.append(fichier)
		return sorted(lsFilesFasta, key=sort_human)

	def splitFilesDir(self):
		"""list files and list directory"""
		self.lsInDir()
		self.listDir = []		# only directory in path
		self.listFiles = []		# only files in path
		# Ouverture des fichiers du repertoire
		for fichier in self.listPath:
			try:
				if os.path.isdir(fichier) == True :		# C'est un dossier
					self.listDir.append(str(fichier))
				elif os.path.exists(fichier) == True :	# C'est un fichier
					self.listFiles.append(str(fichier))
			except:
				print("ERROR MODULES_SEB::Class-directory : path '%s' is not valide path contain other type (not files or directory)" % self.pathDirectory)
				exit()

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
	"""
	if relative[0] != "/":			# The relative path is a relative path, ie do not starts with /
		absolutePath = os.path.realpath(relative)
		return absolutePath
	else:						# Relative is in fact an absolute path, send a warning
		absolutePath = relative;
		return absolutePath

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

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='16/06/2016'
debug="False"
#debug="True"

specFile = """
# -*- mode: python -*-

block_cipher = None

added_files = [ **ADD** ]

a = Analysis(['**PYTHON**'],
             pathex=['**workingDir**'],
             binaries=None,
             datas=added_files,
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='**NAME**',
          debug=False,
          strip=False,
          upx=True,
          console=False , icon='**ICON**', resources=['**UI**'])


"""


##################################################
## Functions

###################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='run_pyinstaller.py', description='''This run pyinstaller to build file on OS''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display run_pyinstaller.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-p', '--python', metavar="<filename>", type=existant_file, required=True, dest = 'pythonFile', help = 'Main programm script')
	filesReq.add_argument('-u', '--ui', metavar="<filename>", type=existant_file, required=True, dest = 'uiFile', help = 'UI file associated to python')
	filesReq.add_argument('-a', '--addIncludes', metavar="<path/to/directory>", type = directory, required=True, dest = 'includePath', help = 'Name of directory with includes file (must be name includes)')
	filesReq.add_argument('-i', '--icon', metavar="<filename>", type=existant_file, required=True, dest = 'iconFile', help = 'icon.ico file')


	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=False, dest = 'pathOut', help = 'Name of output file directory')
	files.add_argument('-py', '--pythonVersion', metavar="2/3", type = int, default=3, required=False, dest = 'pyversion', help = 'python version (default 3)')


	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#              Welcome in run_pyinstaller (Version " + version + ")               #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	## defaults option
	if args.pathOut != None:
		pathOut = args.pathOut

	# Récupère les arguments
	if "win" in sys.platform :
		workingDir = os.path.dirname(os.path.realpath(args.pythonFile))+"\\"
		distPath = workingDir+"dist\\"
		buildPath = workingDir+"build\\"
		txtADD = "('%s\\\*', 'includes')," % args.includePath
		os_run = "windows"
		if args.pathOut == None:
			pathOut = workingDir+"\windows\\"
			print(pathOut)
			if not os.path.exists(pathOut):
				os.mkdir(pathOut)
			else:
				shutil.rmtree(pathOut)
				os.mkdir(pathOut)

	if "linux" in sys.platform :
		workingDir = os.path.dirname(os.path.realpath(args.pythonFile))+"/"
		distPath = workingDir+"dist/"
		buildPath = workingDir+"build/"
		txtADD = "('%s*', 'includes')," % args.includePath
		os_run = "linux"
		if args.pathOut == None:
			pathOut = workingDir+"/linux/"
			if not os.path.exists(pathOut):
				os.mkdir(pathOut)
			else:
				shutil.rmtree(pathOut)
				os.mkdir(pathOut)

	# resume value to user
	printCol.purple(" - Intput Info:")
	printCol.purple("\t - Python main script is: %s" % args.pythonFile)
	printCol.purple("\t - UI file is: %s" % args.uiFile)
	printCol.purple("\t - icon file is: %s" % args.iconFile)
	printCol.purple("\t - all file in directory '%s' were use" % args.includePath)

	printCol.purple(" - Output Info:")
	printCol.purple("\t - Working directory: %s" % workingDir)
	printCol.purple("\t - Output executable were in directory: %s\n\n" % pathOut)

	icon = os.path.basename(args.iconFile)
	ui = os.path.basename(args.uiFile)
	pythonFileBase = os.path.basename(args.pythonFile)
	nameGUI = "%s_%s" % (pythonFileBase.split(".")[0], os_run )

	dictToReplace = {
	"**ADD**"		:	str(txtADD),
	"**ICON**"	:	str(args.iconFile),
	"**UI**"		:	str(ui),
	"**PYTHON**"		:	str(pythonFileBase),
	"**NAME**"		:	str(nameGUI),
	"**workingDir**"		:	str(workingDir.replace("\\","\\\\"))
	}

	specFileReplace = replace_all(dictToReplace, specFile)


	# build file spec
	specFileName = workingDir+nameGUI+".spec"
	with open(specFileName, "w") as specFile:
		specFile.write(specFileReplace)


	# Test if already run to remove old compilation
	if os.path.isdir(distPath):
		shutil.rmtree(distPath)
	if os.path.isdir(buildPath):
		shutil.rmtree(buildPath)

	# run PyInstaller
	if args.pyversion == 2:
		os.system("python -m PyInstaller %s" %specFile.name)
	if args.pyversion == 3:
		os.system("python3 -m PyInstaller %s" %specFile.name)

	# if PyInstaller run correctly
	buildFile = os.listdir(distPath)[0]			# list dist file
	shutil.move(distPath+buildFile, pathOut)	# move to destOut
	shutil.rmtree(distPath)						# remove compilation
	shutil.rmtree(buildPath)					# remove dist

	printCol.purple("\n\nexecutable product on %s" %pathOut+buildFile )




	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
