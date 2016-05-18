#!/usr/bin/python3.4
# -*- coding: utf-8 -*-
## @package renameFile.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
## Python modules
import argparse, os, subprocess, sys
from time import localtime, strftime
import glob, re

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='03/05/2016'

#################################################
# CLASS
#################################################
#*********************************************** Classe directory *******************
class directory(str):
	"""Class which derives from string.
		Checks that the string is and path to valid directory and not empty"

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

	def splitFilesDir(self):
		"""list files and list directory"""
		self.lsInDir()
		self.listDir = []		# only directory in path
		self.listFiles = []		# only files in path
		# Ouverture des fichiers du repertoire
		for fichier in self.listPath:
			try:
				if os.path.isdir(fichier) == True :		# C'est un dossier
					self.listDir.append(fichier)
				elif os.path.exists(fichier) == True :	# C'est un fichier
					self.listFiles.append(fichier)
			except:
				print("ERROR MODULES_SEB::Class-directory : path '%s' is not valide path contain other type (not files or directory)" % self.pathDirectory)
				exit()



##################################################
## Fonctions

def relativeToAbsolutePath(relative):
	"""	Return the absolutPath

	:param relative: a string path
	:type relative: string
	:rtype: string()
	:return: absolutePath
	:warn: need subprocess::check_output

	Example:
		>>> print(relative)
			../test
		>>> pathDirectory = relativeToAbsolutePath(relative)
		>>> print(pathDirectory)
			/home/sebastien/test

	"""
	from subprocess import check_output
	if relative[0] != "/":			# The relative path is a relative path, ie do not starts with /
		command = "readlink -m "+relative
		absolutePath = subprocess.check_output(command, shell=True).decode("utf-8").rstrip()
		return absolutePath
	else:						# Relative is in fact an absolute path, send a warning
		absolutePath = relative;
		return absolutePath

def replace_all(repls, str):
	return re.sub('|'.join(re.escape(key) for key in repls.keys()),lambda k: repls[k.group(0)], str)


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='renameFile.py', description='''This Programme rename files in and directory''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display renameFile version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory', metavar="<path/to/directory>", type = directory, required=True, dest = 'dirPath', help = 'path with files to rename')
	filesreq.add_argument('-r', '--replace', metavar="<OLD:NEW,OLD:NEW>", required=True, dest = 'replaceParam', help = 'Expression replace must be OLD:NEW and if multiple replace use comma to separate ')

	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#         Welcome in renameFile (Version " + version + ")           #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	workingObjDir = args.dirPath
	replaceParam = args.replaceParam

	print("Workink Directory: %s" % workingObjDir.pathDirectory)

	dicoReplace = {}
	for value in replaceParam.split(","):
		old,new = value.split(":")
		dicoReplace[old] = new


	for fileIn in workingObjDir.listFiles:
		basename = fileIn.split("/")[-1].split(".")[0]
		extention = "."+".".join(fileIn.split("/")[-1].split(".")[1:])

		newName = replace_all(dicoReplace, basename)
		print(basename+extention,"\t",newName+extention)
		os.rename(fileIn , workingObjDir.pathDirectory+newName+extention)


	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
