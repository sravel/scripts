#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
# @package run_multiblast.py
# @author Sebastien Ravel

"""
	The run_multiblast script
	=========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 08/07/2016
	:version: 0.1

	Script description
	------------------

	This Programme run 1 blast for file include in path.\n
	If you want split multifasta into one seq/fasta use script splitMultiFasta.py

	Example
	-------

	>>> run_multiblast.py -f fasta/ -of 5 -o XML/

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display run_multiblast.py version number and exit

	Input mandatory infos for running:
		- \-f <path/to/directory/fasta>, --fasta <path/to/directory/fasta>
						path to fasta files
		- \-o <path/to/directory>, --out <path/to/directory>
						Name of output file directory

	Input infos for running with default values:
		- \-t <sting>, --type <sting>
						Type of blast (blastx, blastn, ...) (default = blastx)
		- \-b <path/to/directory/bank>, --databank <path/to/directory/bank>
						Path to bank fasta (default = /work/BANK/biomaj/nr/current/flat/nr)
		- \-of <int/string>, --outfmt <int/string>
						outfmt of blast (default = 6)
		- \-bo [<string> [<string> ...]], --blastoption [<string> [<string> ...]]
						Other blast options like -bo "-evalue 10-3" "-gapopen 5" (default = "")
		- \-j <int>, --nbjob <int>
						Number of job array lunch (default = 100)
		- \-e <string>, --extention <string>
						Extention of blast output file (default = ".txt")

"""

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath, directory, printCol

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
	parser = argparse.ArgumentParser(prog='run_multiblast.py', description='''This Programme run 1 blast for file include in path. if you want split multifasta into one seq/fasta use script splitMultiFasta.py''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display run_multiblast.py version number and exit')
	parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-f', '--fasta', metavar="<path/to/directory/fasta>", type=directory, required=True, dest = 'fastaFileDir', help = 'path to fasta files')
	filesReq.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'pathOut', help = 'Name of output file directory')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-t', '--type', metavar="<sting>",type = str, default="blastx", required=False, dest = 'typeBlast', help = 'Type of blast (blastx, blastn, ...) (default = blastx)')
	files.add_argument('-b', '--databank', metavar="<path/to/directory/bank>", default="/work/BANK/biomaj/nr/current/flat/nr",required=False, dest = 'dbPath', help = 'Path to bank fasta (default = /work/BANK/biomaj/nr/current/flat/nr')
	files.add_argument('-of', '--outfmt', metavar="<int/string>",type = str, default="6",required=False, dest = 'outfmtValue', help = 'outfmt of blast (default = 6)')
	files.add_argument('-bo', '--blastoption', metavar="<string>", nargs='*', default=[""],required=False, dest = 'blastOptionValue', help = 'Other blast options like -bo "-evalue 10-3" "-gapopen 5" (default = "")')
	files.add_argument('-j', '--nbjob', metavar="<int>", type = int, default=100,required=False, dest = 'nbJobValue', help = 'Number of job array lunch (default = 100)')
	files.add_argument('-e', '--extention', metavar="<string>", type = str, default=".txt",required=False, dest = 'extValue', help = 'Extention of blast output file (default = ".txt")')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#           Welcome in run_multiblast (Version " + version + ")            #")
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
	SGENameFile = outputSHDir+"submitQsubBLAST.sge"

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
	print("\t - Output trash were in directory: %s\n\n" % outputTrashDir)


	# build directory out
	os.makedirs(outputSHDir, exist_ok=True)															# création d'un dossier sh_scripts pour lancer les analyses structures
	os.makedirs(outputTrashDir, exist_ok=True)
	os.makedirs(outputBlastResDir, exist_ok=True)

	count = 1
	for fasta in pathFastaFile.lsExtInDirToList(["fasta"]):
		basenameFasta = fasta.split("/")[-1].split(".")[0]

		with open(outputSHDir+str(count)+"_blast.sh", "w") as shScript:
			shScript.write("module load compiler/gcc/4.9.2 bioinfo/ncbi-blast/2.2.30\n")
			blastcmd = "%s -query %s -db %s -outfmt %s %s -out %s" % (typeBlast, fasta, dbPath, outfmtValue, blastOptionValue, outputBlastResDir+basenameFasta+args.extValue)
			if args.debug == "True" : print(blastcmd)
			shScript.write(blastcmd)
		count+=1



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

/bin/bash """+outputSHDir+"""${SGE_TASK_ID}_blast.sh"""


	with open(SGENameFile, "w") as SGEFile:
		SGEFile.write(headerSGE)

	print("\n - Execution summary:")

	print("\n  You want run MutilblastX for %s fasta,\
 The script are created all fasta-MutilblastX.sh for all fasta into %s,\n\
 For run all sub-script in qsub, %s was created.\n\
 It lunch programm with job array and run %s job max:\n" %(count-1,outputSHDir,SGENameFile, args.nbJobValue))
	printCol.green("\tqsub %s" % SGENameFile)
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
