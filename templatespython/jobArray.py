#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package scriptName.py
# @author Sebastien Ravel

"""
	The scriptName script
	=========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 08/07/2016
	:version: 0.1

	Script description
	------------------

	This Programme ...

	Example
	-------

	>>> scriptName.py -f fasta/ -o resultDir -j 100

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display scriptName.py version number and exit

	Input mandatory infos for running:
		- \-f <path/to/directory/fasta>, --fasta <path/to/directory/fasta>
						path to fasta files
		- \-o <path/to/directory>, --out <path/to/directory>
						Name of output file directory

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
from MODULES_SEB import directory, printCol

## Python modules
import argparse
from time import localtime, strftime

## BIO Python modules


##################################################
## Variables Globales
version="0.1"
VERSION_DATE='27/01/2017'


##################################################
## Functions

###################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme ....''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug', action='store_true', dest='debug', help='enter verbose/debug mode')

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-f', '--fasta', metavar="<path/to/directory>", type=directory, required=True, dest = 'filesDir', help = 'Path to loop files')
	filesReq.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'pathOut', help = 'Name of output file directory (must exist)')

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
	pathFiles = args.filesDir
	pathFilesOut = args.pathOut

	# build directory out
	outputResDir = pathFilesOut.pathDirectory+"scriptRes/"
	outputSHDir = pathFilesOut.pathDirectory+"sh/"
	outputTrashDir = pathFilesOut.pathDirectory+"trash/"
	SGENameFile = outputSHDir+"submitQsubBLAST.sge"

	os.makedirs(outputSHDir, exist_ok=True)
	os.makedirs(outputTrashDir, exist_ok=True)
	os.makedirs(outputResDir, exist_ok=True)

	# resume value to user
	print(" - Intput Info:")
	print("\t - Working in directory: %s" % pathFilesOut.pathDirectory)
	print("\t - Fasta were in directory: %s" % pathFiles.pathDirectory)

	print(" - Output Info:")
	print("\t - Output with result were in directory: %s" % outputResDir)
	print("\t - Output sh were in directory: %s" % outputSHDir)
	print("\t - Output trash were in directory: %s\n\n" % outputTrashDir)



	count = 1
	for fastaFile in pathFiles.lsExtInDirToList(["fasta","fa","fn","fna"]): 					# change with you own extention file
		basenameFasta = fastaFile.split("/")[-1].split(".")[0]									# if
		if args.debug : print(basenameFasta)

		with open(outputSHDir+str(count)+"_script.sh", "w") as shScript:						# not change value here

			shScript.write("module load compiler/gcc/4.9.2 bioinfo/ncbi-blast/2.2.30 \n")		# add your own module load
			cmd = "yourProgramme -f %s -o %s \n" % (fastaFile, outputResDir+basenameFasta)	# add your programme command line
			shScript.write(cmd)

			if args.debug : print(cmd)
		count+=1

	# Not change anythings for headerSGE
	headerSGE = """
#!/bin/bash

#$ -N blast
#$ -cwd
#$ -V
#$ -e """+outputTrashDir+"""
#$ -o """+outputTrashDir+"""
#$ -q long.q
#$ -t 1-"""+str(count-1)+"""
#$ -tc """+str(args.nbJobValue)+"""
#$ -S /bin/bash

/bin/bash """+outputSHDir+"""${SGE_TASK_ID}_script.sh"""

	# write the .sge file
	with open(SGENameFile, "w") as SGEFile:
		SGEFile.write(headerSGE)


	# print summary, you can modifie but keep printCol line
	print("\n - Execution summary:")
	print("\n You want run script for %s files,\
 The script are created all ${SGE_TASK_ID}_script.sh for all files into %s,\n\
 For run all sub-script in qsub, %s was created.\n\
 It lunch programm with job array and run %s job max:\n" %(count-1,outputSHDir,SGENameFile, args.nbJobValue))
	printCol.green("\tqsub %s" % SGENameFile)


	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
