#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
# @package compress.py
# @author Sebastien Ravel

"""
	The compress script
	=========================

	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 21/10/2016
	:version: 0.1

	Script description
	------------------

	This Programme compress all file with pass extention in directory with job array

	Example
	-------

	>>> compress.py -f ./fasta/ -e fastq

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display compress.py version number and exit

	Input mandatory infos for running:
		- \-f <path/to/directory/fasta>, --fasta <path/to/directory/fasta>
						path to fasta files
		- \-e string [string ...], --extention string [string ...]
						list of extention to compress

	Input infos for running with default values:
		- \-j <int>, --nbjob <int>
						Number of job array lunch (default = 100)
"""

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
sys.path.insert(1,'/homedir/sravel/programme/ScriptsSEB/scripts//modules/')
from MODULES_SEB import relativeToAbsolutePath, directory, printCol

## Python modules
import argparse
from time import localtime, strftime

## BIO Python modules


##################################################
## Variables Globales
version="0.1"
VERSION_DATE='21/10/2016'
debug="False"
#debug="True"


##################################################
## Functions

###################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme compress all file with pass extention in directory with job array''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-f', '--fasta', metavar="<path/to/directory/fasta>", type=directory, required=True, dest = 'fastaFileDir', help = 'path to fasta files')
	#filesReq.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'pathOut', help = 'Name of output file directory')
	filesReq.add_argument('-e', '--extention', metavar="string", type = str, nargs='+', required=True, dest = 'extentionList', help = 'list of extention to compress')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-j', '--nbjob', metavar="<int>", type = int, default=100,required=False, dest = 'nbJobValue', help = 'Number of job array lunch (default = 100)')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#           Welcome in %s (Version %s)            #" %(__file__, version))
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les arguments
	pathFileFile = args.fastaFileDir
	#pathFileOut = args.pathOut
	pathFileOut = pathFileFile


	outputSHDir = pathFileOut.pathDirectory+"sh/"
	outputTrashDir = pathFileOut.pathDirectory+"trash/"
	SGENameFile = outputSHDir+"submitQsubCompress.sge"


	# resume value to user
	print(" - Intput Info:")
	print("\t - Files were in directory: %s" % pathFileFile.pathDirectory)
	print("\t - Extention to compress: %s" % args.extentionList)


	print(" - Output Info:")
	print("\t - Output directory: %s" % pathFileOut.pathDirectory)
	print("\t - Output sh were in directory: %s" % outputSHDir)
	print("\t - Output trash were in directory: %s\n\n" % outputTrashDir)


	# build directory out
	os.makedirs(outputSHDir, exist_ok=True)															# création d'un dossier sh_scripts pour lancer les analyses structures
	os.makedirs(outputTrashDir, exist_ok=True)


	count = 1
	for fasta in pathFileFile.lsExtInDirToList(args.extentionList):
		print("FASTAFILE:",fasta)
		basenameFasta = fasta.split("/")[-1]
		basedir = "/".join(fasta.split("/")[:-1])

		with open(outputSHDir+str(count)+"_tar.sh", "w") as shScript:

			#shScript.write("""# Creates an archive (*.tar.gz) from given directory.\nfunction maketar() { tar cvzf "${1%%/}.tar.gz"  "${1%%/}/"; }\n\n""")
			#shScript.write("cd %s\nmaketar %s\n" % (basedir,basenameFasta))
			shScript.write("cd %s\ngzip %s\n" % (basedir,basenameFasta))
		count+=1



	headerSGE = """
#!/bin/bash

#$ -N compress
#$ -cwd
#$ -V
#$ -e """+outputTrashDir+"""
#$ -o """+outputTrashDir+"""
#$ -q long.q
#$ -t 1-"""+str(count-1)+"""
#$ -tc """+str(args.nbJobValue)+"""
#$ -S /bin/bash

/bin/bash """+outputSHDir+"""${SGE_TASK_ID}_tar.sh"""


	with open(SGENameFile, "w") as SGEFile:
		SGEFile.write(headerSGE)

	print("\n - Execution summary:")

	print("\n  You want compress %s files,\
 The script are created all {}_tar.sh for all files into %s,\n\
 For run all sub-script in qsub, %s was created.\n\
 It lunch programm with job array and run %s job max:\n" %(count-1,outputSHDir,SGENameFile, args.nbJobValue))
	printCol.green("\tqsub %s" % SGENameFile)
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
