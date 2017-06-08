#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package buildJobArray.py
# @author Sebastien Ravel

"""
	The buildJobArray script
	=========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 08/06/2017
	:version: 0.1

	Script description
	------------------

	This Programme run bash script in buildJobArray

		- \-directory 'fasta' contains x fasta files with extention fasta.

		- \-script.sh is:

	>>> #!/bin/bash
	>>> module load compiler/gcc/4.9.2 bioinfo/ncbi-blast/2.2.30
	>>> blastx -query FILEIN -db /work/BANK/biomaj/nr/current/flat/nr -outfmt 5  -out BASENAMEOUTEXTOUT

	By default script buildJobArray.py replace
		- FILEIN by file name in loop
		- BASENAMEOUT by basePath+basename (basename is name of file in loop (use only one . in file name)

	Example
	-------

	>>> buildJobArray.py -f fasta/ -e fasta -o ./ -s script.sh -r EXTOUT,.blast

	Help Programm
	-------------

	optional arguments:
		- \-h, --help            show this help message and exit
		- \-v, --version         display buildJobArray.py version number and exit
		- \-dd, --debug          enter verbose/debug mode

	Input mandatory infos for running:
		- \-f <path/to/directory>, --files <path/to/directory>
				Path to loop files
		- \-e str, --extention str
				extention for select files to loop
		- \-o <path/to/directory>, --out <path/to/directory>
				Name of output file directory (must exist)
		- \-s <path/to/shScriptOut>, --sh <path/to/shScriptOut>
				Name sh file with programme loop for jobArray
		- \-r FILEIN:fileName [FILEIN:fileName ...], --replace FILEIN:fileName [FILEIN:fileName ...]
				Replace any word in sh script (from -s) by second
				value, for example: EXTENTIONOUT,.txt

	nput infos for running with default values:
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
from MODULES_SEB import directory, printCol, replace_all, existant_file

## Python modules
import argparse
from time import localtime, strftime

def listTuple(s):
	try:
		x, y = map(str, s.split(','))
		return x, y
	except Exception:
		raise argparse.ArgumentTypeError("'{0}' does not with correct format, Must be 'value1,value2'.".format(s))


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
	filesReq.add_argument('-f', '--files', metavar="<path/to/directory>", type=directory, required=True, dest = 'filesDir', help = 'Path to loop files')
	filesReq.add_argument('-e', '--extention', metavar="str", type=str, required=True, dest = 'fileExt', help = 'extention for select files to loop')
	filesReq.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'pathOut', help = 'Name of output file directory (must exist)')
	filesReq.add_argument('-s', '--sh', metavar="<path/to/shScriptOut>", type = existant_file, required=True, dest = 'shFile', help = 'Name sh file with programme loop for jobArray')
	filesReq.add_argument('-r', '--replace', metavar="FILEIN:fileName", type=listTuple, dest="listTuple", nargs="+", help="Replace any word in sh script (from -s) by second value, for example: EXTENTIONOUT,.txt")

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
	shFile = args.shFile
	fileExt = args.fileExt
	dictToReplace = dict(args.listTuple)

	# build directory out
	outputResDir = pathFilesOut.pathDirectory+"scriptRes/"
	outputSHDir = pathFilesOut.pathDirectory+"sh/"
	outputTrashDir = pathFilesOut.pathDirectory+"trash/"
	SGENameFile = outputSHDir+"submitQsubJOBARRAY.sge"

	os.makedirs(outputSHDir, exist_ok=True)
	os.makedirs(outputTrashDir, exist_ok=True)
	os.makedirs(outputResDir, exist_ok=True)

	# resume value to user
	print(" - Intput Info:")
	print("\t - Working in directory: %s" % pathFilesOut.pathDirectory)
	print("\t - Files were in directory: %s" % pathFiles.pathDirectory)
	print("\t - FileExt is is  %s" % fileExt)
	print("\t - SH file is  %s" % shFile)

	print(" - Output Info:")
	print("\t - Output with result were in directory: %s" % outputResDir)
	print("\t - Output sh were in directory: %s" % outputSHDir)
	print("\t - Output trash were in directory: %s\n\n" % outputTrashDir)

	with open(shFile,"r") as shFileIn:
		shtxt = "".join(shFileIn.readlines())


	count = 1
	for fileName in pathFiles.lsExtInDirToList([fileExt]): 					# change with you own extention file
		basenameFile = fileName.split("/")[-1].split(".")[0]


		if args.debug : print("############ loop "+str(count)+" ############")
		if args.debug : print(">>> FILEIN: "+fileName)
		if args.debug : print(">>> BASENAMEOUT: "+outputResDir+"/"+basenameFile)

		dictToReplace.update({
			"FILEIN"		:	str(fileName),
			"BASENAMEOUT"	:	str(outputResDir+"/"+basenameFile),
			})

		with open(outputSHDir+str(count)+"_script.sh", "w") as shScriptOut:						# not change value here

			shtxtReplace = replace_all(dictToReplace, shtxt )
			shScriptOut.write(shtxtReplace)

			if args.debug : print(">>> base Script with modifications: \n")
			if args.debug : print(shtxtReplace)
		count+=1

	# Not change anythings for headerSGE
	headerSGE = """
#!/bin/bash

#$ -N SEB_buildJobArray
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
	print("\n ######################################")
	print("\n - Execution summary:")
	print("\n You want run script for %s files,\
 The script are created all ${SGE_TASK_ID}_script.sh for all files into %s,\n\
 For run all sub-script in qsub, %s was created.\n\
 It run programm with job array with %s job max:\n" %(count-1,outputSHDir,SGENameFile, args.nbJobValue))
	printCol.green("\tqsub %s" % SGENameFile)


	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
