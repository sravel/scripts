#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package buildTracks.py
# @author Sebastien Ravel

"""
	The buildTracks script
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

	>>> buildTracks.py -f fasta/ -r referance.fasta -o output/

	Help Programm
	-------------

	optional arguments:
		- \-h, --help            show this help message and exit
		- \-v, --version         display ./buildTracks.py version number and exit
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
VERSION_DATE='08/11/2016'


###################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme build track file for jBrowse''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug',action='store_true', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-d', '--dir', metavar="<path/to/directory/>", type=directory, required=True, dest = 'fastaFileDir', help = 'path to files files')
	filesReq.add_argument('-e', '--extention', metavar="<string>", choices=["bam", "vcf"], required=True, dest = 'extentionFile', help = 'file you want to build')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-k', '--key', metavar="<string>", type = str, default="",required=False, dest = 'keyInfo', help = 'name of track')
	files.add_argument('-c', '--category', metavar="<string>",type = str, default="Mapping", required=False, dest = 'catInfo', help = 'name of category of tracks')
	files.add_argument('-o', '--output', metavar="<path/to/file>",type = existant_file, default="trackFile", required=False, dest = 'trackFile', help = 'name of output file with tracks')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#           Welcome in %s (Version %s)            #" %(__file__, version))
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les arguments
	pathFastaFile = args.fastaFileDir
	extentionFile = args.extentionFile

	# defaults option
	keyInfo = args.keyInfo
	catInfo = args.catInfo
	trackFileName = args.trackFile

	# resume value to user
	print(" - Intput Info:")
	print("\t - Files are in directory: %s" % pathFastaFile.pathDirectory)
	print("\t - ExtentionFile is : %s" % extentionFile)

	print(" - Output Info:")
	print("\t - Output file is: %s" % trackFileName)
	print("\t - catInfo file is: %s" % catInfo)
	print("\t - keyInfo file is: %s" % keyInfo)



	pathFastaFile = directory("/work/carlier.j/assemblyFastq0615/")

	with open(trackFileName,"w") as trackFile:
		for fasta in pathFastaFile.lsExtInDirToList(extentionFile):
			basenameDir = fasta.split("/")[-1].split(".")[0]
			print(basenameDir)

			if extentionFile == "bam":
				strTrack = """
[ tracks.{0}]
storeClass     = JBrowse/Store/SeqFeature/BAM
urlTemplate    = {1}/{0}.{1}
baiUrlTemplate = {1}/{0}.bai
category = {4}
type = JBrowse/View/Track/Alignments2
key  = {3}{0}
""".format(basenameDir, extentionFile, keyInfo, catInfo)

				trackFile.write(strTrack)

			if extentionFile == "vcf":
				strTrack = """
[ tracks.{0}]
storeClass     = JBrowse/Store/SeqFeature/VCFTabix
urlTemplate    = {1}/{0}.vcf.gz
tbiUrlTemplate = {1}/{0}.vcf.gz.tbi
category = {4}
type = JBrowse/View/Track/CanvasVariants
key  = {3}{0}
""".format(basenameDir, extentionFile, keyInfo, catInfo)

				trackFile.write(strTrack)



















	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
