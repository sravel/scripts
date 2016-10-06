#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
# @package run_Assembly.py
# @author Sebastien Ravel
#__docformat__ = "restructuredtext en"

"""
	The run_Assembly script
	=======================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 03/05/2016
	:version: 0.1

	Script description
	------------------

	This Programme run assembly of Jerome Gouzi pipeline for fastq file

	Example
	-------

	>>> run_Assembly.py -f fastq/ -o assemblyFinal

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display run_Assembly version number and exit

	Input mandatory infos for running:
		- \-f <path/to/directory>, --directory <path/to/directory>
						path with Fastq files
		- \-o <path/to/directory>, --out <path/to/directory>
						Output Path

"""

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import directory

## Python modules
import argparse
from time import localtime, strftime


##################################################
## Variables Globales
version="0.1"
VERSION_DATE='03/05/2016'


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='run_Assembly.py', description='''This Programme run assembly of Jerome Gouzi pipeline for fastq file, fastq MUST be ?z format''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display run_Assembly version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fastq', metavar="<path/to/directory>", type = directory, required=True, dest = 'dirPath', help = 'path With Fastq files')
	filesreq.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'outPath', help = 'Output Path')

	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#            Welcome in run_Assembly (Version " + version + ")              #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	pathFastqFile = args.dirPath
	pathFileOut = args.outPath

	outputSHDir = pathFileOut.pathDirectory+"sh/"
	outputTrashDir = pathFileOut.pathDirectory+"trash/"
	SGENameFile = outputSHDir+"submitQsubBLAST.sge"


	# resume value to user
	print(" - Intput Info:")
	print("\t - Working in directory: %s" % pathFileOut.pathDirectory)
	print("\t - Fastq were in directory: %s" % pathFastqFile.pathDirectory)

	print(" - Output Info:")
	print("\t - Output were in directory: %s" % pathFileOut.pathDirectory)
	print("\t - Output sh were in directory: %s" % outputSHDir)
	print("\t - Output trash were in directory: %s\n\n" % outputTrashDir)


	# build directory out
	os.makedirs(outputSHDir, exist_ok=True)															# création d'un dossier sh_scripts pour lancer les analyses structures
	os.makedirs(outputTrashDir, exist_ok=True)

	nbJob = str(int(len(pathFastqFile.listFiles)/2))

	# Parcours des fichiers
	count=1
	listFiles = []

	for fileIn in pathFastqFile.listFiles:
		cmd = ""
		txt = ""
		fileName = fileIn.split("/")[-1]
		extention = fileName.split(".")[-1]
		baseDir = "/".join(fileIn.split("/")[:-1])
		#basename = fileIn.split("/")[-1].split(".")[0].split("_")[0]
		basename = fileIn.split("/")[-1].split(".")[0]
		rValue = "_"+fileIn.split("/")[-1].split(".")[0].split("_")[-1]
		extention = "."+".".join(fileIn.split("/")[-1].split(".")[1:])
		print(basename+rValue+extention)

		if "z" not in extention:
			cmd = """function maketar() { tar cvzf "${1%%/}.tar.gz"  "${1%%/}/"; }\n\n"""
			cmd += "maketar %s*\n\n" % pathFileOut.pathDirectory+basename+"/"+fileName[:-1]
			#print("Warning, file %s must be fastq.?z format" % fileIn)
			#exit()

		os.mkdir(pathFileOut.pathDirectory+basename)

		cmd = "ln -s "+pathFastqFile.pathDirectory+fileName+" "+pathFileOut.pathDirectory+basename+"/"+fileName
		os.system(cmd)

		if basename not in listFiles:
			listFiles.append(basename)

			txt += cmd
			txt += """/NAS/BAILLARGUET/BGPI/tools/lipm_assembly/bin/lipm_assemble_solexa_pe.pl --datadir """+pathFileOut.pathDirectory+basename+""" --outdir """+pathFileOut.pathDirectory+basename+""" --outprefix """+basename+""" --log """+pathFileOut.pathDirectory+basename+"""/log.txt\n"""
			with open(outputSHDir+"/"+str(count)+"-assembly.sh","w") as shScript:
				shScript.write(txt)
			count+=1

	# création d'un script de qsub array
	headerJoB = """
#!/bin/bash

#$ -N ASSEMBLY
#$ -cwd
#$ -V
#$ -e """+outputTrashDir+"""
#$ -o """+outputTrashDir+"""
#$ -q long.q
#$ -t 1-"""+nbJob+"""
#$ -tc 10
#$ -S /bin/bash

qsub -N denovo -q long.q -e trash -o trash -l mem_free=100G -cwd """+outputSHDir+"""/$SGE_TASK_ID-assembly.sh
	"""
	with open(pathFileOut.pathDirectory+"submitQsub.sge","w") as SGEFile:
		SGEFile.write(headerJoB)

	print("\n  - Outputting \n\
\t- %i directories have been created base on the file name, contain link to fastq and output of assembly\n\
\t- Two directories were created : 'sh' and 'trash', sh contain commande lines and trash all the .o of jobs\n\
\n  - To launch all assembly execute on the cluster: \n\
\tqsub %s\n" %(count-1,SGEFile.name))

	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
