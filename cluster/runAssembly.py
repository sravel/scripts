#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
## @package runAssembly.py
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


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='runAssembly.py', description='''This Programme rename files in and directory''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display runAssembly version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory', metavar="<path/to/directory>", type = directory, required=True, dest = 'dirPath', help = 'path With Fastq files')
	filesreq.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'outPath', help = 'Output Path')

	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#             Welcome in runAssembly (Version " + version + ")              #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	workingObjDir = args.dirPath
	outObjDir = args.outPath

	print("Workink Directory: %s" % workingObjDir.pathDirectory)

	# création d'un output assembly dans le out et trash
	try:
		os.mkdir(outObjDir.pathDirectory+"trash")
	except FileExistsError:
		pass
	try:
		os.mkdir(outObjDir.pathDirectory+"sh")
	except FileExistsError:
		pass


	nbJob = str(int(len(workingObjDir.listFiles)/2))

	# création d'un script de qsub array
	headerJoB = """
#!/bin/bash

#$ -N ASSEMBLY
#$ -cwd
#$ -V
#$ -e """+outObjDir.pathDirectory+"trash"+"""
#$ -o """+outObjDir.pathDirectory+"trash"+"""
#$ -q long.q
#$ -t 1-"""+nbJob+"""
#$ -tc 5
#$ -S /bin/bash

qsub -N denovo -q long.q -e trash -o trash -l mem_free=100G -cwd """+outObjDir.pathDirectory+"sh"+"""/$SGE_TASK_ID-assembly.sh
	"""
	SGEFile = open(outObjDir.pathDirectory+"submitQsub.sge","w")
	SGEFile.write(headerJoB)
	#print(headerJoB)

	# Parcours des fichiers
	count=1
	listFiles = []
	for fileIn in workingObjDir.listFiles:
		fileName = fileIn.split("/")[-1]
		basename = fileIn.split("/")[-1].split(".")[0].split("_")[0]
		rValue = "_"+fileIn.split("/")[-1].split(".")[0].split("_")[-1]
		extention = "."+".".join(fileIn.split("/")[-1].split(".")[1:])

		#print(basename+rValue+extention)
		try:
			os.mkdir(outObjDir.pathDirectory+basename)
		except FileExistsError:
			pass
		cmd = "ln -s "+workingObjDir.pathDirectory+fileName+" "+outObjDir.pathDirectory+basename+"/"+fileName
		os.system(cmd)

		if basename not in listFiles:
			listFiles.append(basename)

			txt = """/NAS/BAILLARGUET/BGPI/tools/lipm_assembly/bin/lipm_assemble_solexa_pe.pl --datadir """+outObjDir.pathDirectory+basename+""" --outdir """+outObjDir.pathDirectory+basename+""" --outprefix """+basename+""" --log """+outObjDir.pathDirectory+basename+"""/log.txt --mink 50 --minlen 200\n"""
			with open(outObjDir.pathDirectory+"sh"+"/"+str(count)+"-assembly.sh","w") as shScript:
				shScript.write(txt)
			count+=1
	SGEFile.close()

	print("\n  - Outputting \n\
\t- %i directories have been created base on the file name, contain link to fastq and output of assembly\n\
\t- Two directories were created : 'sh' and 'trash', sh contain commande lines and trash all the .o of jobs\n\
\n  - To launch all assembly execute on the cluster: \n\
\tqsub %s\n" %(count-1,SGEFile.name))

	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
