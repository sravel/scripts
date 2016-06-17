#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
## @package run_multiblastX.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath, existant_file, directory, printCol

## Python modules
import argparse
from time import localtime, strftime

## BIO Python modules


##################################################
## Variables Globales
version="0.1"
VERSION_DATE='16/06/2016'
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
	parser = argparse.ArgumentParser(prog='run_multiblastX.py', description='''This Programme extract Fasta Seq with liste keep''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display run_multiblastX.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-f', '--fasta', metavar="<path/to/directory/fasta>", type=directory, required=True, dest = 'fastaFileDir', help = 'path to fasta files')
	filesReq.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'pathOut', help = 'Name of output file directory')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-t', '--type', metavar="<sting>",type = str, default="blastx", required=False, dest = 'typeBlast', help = 'Type of blast (blastx, blastn, ...) (default = blastx)')
	files.add_argument('-b', '--databank', metavar="<path/to/directory/bank>", default="/work/BANK/biomaj/nr/nr_2016-05-21/flat/nr",required=False, dest = 'dbPath', help = 'Path to bank fasta (default = /work/BANK/biomaj/nr/nr_2016-05-21/flat/nr)')
	files.add_argument('-of', '--outfmt', metavar="<int/string>",type = str, default="6",required=False, dest = 'outfmtValue', help = 'outfmt of blast (default = 6)')
	files.add_argument('-bo', '--blastoption', metavar="<string>", nargs='*', default=[""],required=False, dest = 'blastOptionValue', help = 'Other blast options like -bo "-evalue 10-3" "-gapopen 5" (default = "")')


	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#        Welcome in run_multiblastX (Version " + version + ")          #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les arguments
	pathFastaFile = args.fastaFileDir
	pathFileOut = args.pathOut

	# defaults option
	typeBlast=args.typeBlast
	dbPath=relativeToAbsolutePath(args.dbPath)
	outfmtValue=args.outfmtValue
	blastOptionValue=" ".join(args.blastOptionValue)

	outputBlastResDir = pathFileOut.pathDirectory+"blastRes/"
	outputSHDir = pathFileOut.pathDirectory+"sh/"
	outputTrashDir = pathFileOut.pathDirectory+"trash/"
	SGENameFile = outputSHDir+"submitQsub.sge"

	# resume value to user
	print(" - Intput Info:")
	print("\t - Working in directory: %s" % pathFileOut.pathDirectory)
	print("\t - Fasta were in directory: %s" % pathFastaFile.pathDirectory)
	print("\t - blats is: %s" % typeBlast)
	print("\t - dataBase is: %s" % dbPath)
	print("\t - outfmtValue is: %s" % outfmtValue)
	print("\t - Other options are: %s" % blastOptionValue)

	print(" - Output Info:")
	print("\t - Output with result Blast were in directory: %s" % outputBlastResDir)
	print("\t - Output sh were in directory: %s" % outputSHDir)
	print("\t - Output trash were in directory: %s" % outputTrashDir)


	# build directory out
	os.makedirs(outputSHDir, exist_ok=True)															# création d'un dossier sh_scripts pour lancer les analyses structures
	os.makedirs(outputTrashDir, exist_ok=True)
	os.makedirs(outputBlastResDir, exist_ok=True)

	count = 1
	for fasta in pathFastaFile.lsExtInDirToList("fasta"):
		basenameFasta = fasta.split("/")[-1].split(".")[0]

		with open(outputSHDir+count+"_blast.sh", "w") as shScript:
			shScript.write("module load compiler/gcc/4.9.2 bioinfo/ncbi-blast/2.2.30")
			blastcmd = "%s -query %s -db %s -outfmt %s %s -out %s" % (typeBlast, fasta, dbPath, outfmtValue, blastOptionValue, outputBlastResDir+basenameFasta+".txt")
			print(blastcmd)
			shScript.write(blastcmd)



	headerSGE = """
#!/bin/bash

	#$ -N blast
	#$ -cwd
	#$ -V
	#$ -e """+outputTrashDir+"""
	#$ -o """+outputTrashDir+"""
	#$ -q long.q
	#$ -t 1-"""+count+"""
	#$ -tc 100
	#$ -S /bin/bash

	/bin/bash """+outputSHDir+"""$SGE_TASK_ID'_blast.sh"""


	with open(SGENameFile, "w") as SGEFile:
		SGEFile.write(headerSGE)

	print("\n\nExecution summary:")

	print("\n\n You want run MutilblastX for %s fasta,\
 The script are created all fasta-MutilblastX.sh for all fasta into %s,\n\
 For run all sub-script in qsub, %s was created, It lunch programm make:\n\n" %(count,outputSHDir,SGENameFile))
	printCol.GREEN("\tqsub %s" % SGENameFile)
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
