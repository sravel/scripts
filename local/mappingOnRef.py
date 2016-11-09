#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
# @package mappingOnRef.py
# @author Sebastien Ravel

"""
	The mappingOnRef script
	=========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 08/11/2016
	:version: 0.1

	Script description
	------------------

	This Programme map fasta with bwa mem on reference file

	Example
	-------

	>>> mappingOnRef.py -f fasta/ -r referance.fasta -o output/

	Help Programm
	-------------

	optional arguments:
		- \-h, --help            show this help message and exit
		- \-v, --version         display ./mappingOnRef.py version number and exit
		- \-dd, --debug          enter verbose/debug mode

	Input mandatory infos for running:
		- \-f <path/to/directory/fasta>, --fasta <path/to/directory/fasta>
						path to fasta files
		- \-r <path/to/dreference>, --ref <path/to/dreference>
						path to reference fasta files
		- \-o <path/to/directory>, --out <path/to/directory>
						Name of output file directory

	Input infos for running with default values:
		- \-j <int>, --nbjob <int>
						Number of job array lunch (default = 100)
		- \-th <int>, --thread <int>
						number of threads for mapping (default = 4)
"""

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath, directory, printCol, existant_file

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
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme map fasta with bwa mem on reference file''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug',action='store_true', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-f', '--fasta', metavar="<path/to/directory/fasta>", type=directory, required=True, dest = 'fastaFileDir', help = 'path to fasta files')
	filesReq.add_argument('-r', '--ref', metavar="<path/to/dreference>", type=existant_file, required=True, dest = 'refFile', help = 'path to reference fasta files')
	filesReq.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'pathOut', help = 'Name of output file directory')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-j', '--nbjob', metavar="<int>", type = int, default=100,required=False, dest = 'nbJobValue', help = 'Number of job array lunch (default = 100)')
	files.add_argument('-th', '--thread', metavar="<int>",type = int, default=4, required=False, dest = 'nbThreads', help = 'number of threads for mapping (default = 4)')


	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#           Welcome in %s (Version %s)            #" %(__file__, version))
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les arguments
	pathFastaFile = args.fastaFileDir
	pathFileOut = args.pathOut
	refFile = args.refFile

	# defaults option
	nbThreads=args.nbThreads

	outputMappingResDir = pathFileOut.pathDirectory+"mappingRes/"
	outputSHDir = pathFileOut.pathDirectory+"sh/"
	outputTrashDir = pathFileOut.pathDirectory+"trash/"
	SGENameFile = outputSHDir+"submitQsubBLAST.sge"


	# resume value to user
	print(" - Intput Info:")
	print("\t - Working in directory: %s" % pathFileOut.pathDirectory)
	print("\t - Fasta are in directory: %s" % pathFastaFile.pathDirectory)
	print("\t - Reference is : %s" % refFile)
	print("\t - Number of threads are: %s" % nbThreads)

	print(" - Output Info:")
	print("\t - Output with result Mapping were in directory: %s" % outputMappingResDir)
	print("\t - Output sh were in directory: %s" % outputSHDir)
	print("\t - Output trash were in directory: %s\n\n" % outputTrashDir)


	# build directory out
	os.makedirs(outputSHDir, exist_ok=True)															# création d'un dossier sh_scripts pour lancer les analyses structures
	os.makedirs(outputTrashDir, exist_ok=True)
	os.makedirs(outputMappingResDir, exist_ok=True)

	count = 1
	for fasta in pathFastaFile.lsExtInDirToList(["fasta","fa","fn","fna"]):
		basenameFasta = fasta.split("/")[-1].split(".")[0]

		with open(outputSHDir+str(count)+"_mapping.sh", "w") as shScript:
			mappingcmd = """
module load bioinfo/bwa/0.7.15 system/java/jdk8 bioinfo/samtools/1.3
"""
			#mappingcmd = """
#bwa aln -n 5  -f %s.sai  %s %s
#""" % (outputMappingResDir+basenameFasta,refFile, fasta)

			#mappingcmd += """
#bwa samse -f %s.sam %s %s.sai %s
#""" % (outputMappingResDir+basenameFasta, refFile, outputMappingResDir+basenameFasta, fasta)

			mappingcmd += """
bwa mem %s %s > %s.sam
""" % (refFile, fasta, outputMappingResDir+basenameFasta)


			mappingcmd += """
java -Xmx8G -jar /usr/local/bioinfo/picard-tools/1.130/picard.jar SortSam VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate CREATE_INDEX=TRUE  INPUT=%s.sam OUTPUT=%s.bam
""" % (outputMappingResDir+basenameFasta, outputMappingResDir+basenameFasta)

			mappingcmd += """
samtools flagstat {0}.bam > {0}.samtoolsFlagstat

rm {0}.sai
rm {0}.sam

""".format(outputMappingResDir+basenameFasta)

			if args.debug == True : print(mappingcmd)

			shScript.write(mappingcmd)
		count+=1


	headerSGE = """
#!/bin/bash

#$ -N mapping
#$ -cwd
#$ -V
#$ -e """+outputTrashDir+"""
#$ -o """+outputTrashDir+"""
#$ -q long.q
#$ -pe parallel_smp """+str(nbThreads)+"""
#$ -t 1-"""+str(count-1)+"""
#$ -tc """+str(args.nbJobValue)+"""
#$ -S /bin/bash

/bin/bash """+outputSHDir+"""${SGE_TASK_ID}_mapping.sh"""


	with open(SGENameFile, "w") as SGEFile:
		SGEFile.write(headerSGE)

	print("\n - Execution summary:")

	print("\n  You want run MutilmappingX for %s fasta,\
 The script are created all fasta-MutilmappingX.sh for all fasta into %s,\n\
 For run all sub-script in qsub, %s was created.\n\
 It lunch programm with job array and run %s job max:\n" %(count-1,outputSHDir,SGENameFile, args.nbJobValue))
	printCol.green("\tqsub %s" % SGENameFile)
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
