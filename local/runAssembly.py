#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
## @package runAssembly.py
# @author Sebastien Ravel
#__docformat__ = "restructuredtext en"
"""
	The runAssembly script
	======================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 18/05/2015
	:version: 2.0

	Use it to transform tabular matrice of SNP on fasta file

	Example:

	>>> python runAssembly.py -t

	Input require
	-------------

	- Name of tab file

	Required Module install
	-----------------------

	This module run with Python 3.x and not Python 2.x

	** Include module of python:
		- argparse, time

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
