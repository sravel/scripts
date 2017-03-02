#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
## @package filterTabToSNP.py
# @author Sebastien Ravel
"""
	The compress script
	=========================

	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 28/02/2017
	:version: 0.1

	Script description
	------------------

	Programme extract U from Tab file and build 3 tab:
	- prefilter: without N for all samples
	- withoutN: not missing data allow for all samples
	- withoutNandR: only line with one or more SNP without N


	Example
	-------

	>>> /filterTabToSNP.py -t ./ -o ./

	Help Programm
	-------------

	usage: ./filterTabToSNP.py [-h] [-v] -t <path/to/tabFileDir>
						[-o <path/to/outputDir>]

	This Programme extract U from Tab file

	optional arguments:
		- \-h, --help
					show this help message and exit
		- \-v, --version
					display ./filterTabToSNP.py version number and exit

	Input mandatory infos for running:
		- \-t <path/to/tabFileDir>, --tab <path/to/tabFileDir>
					path to file tab (gzip or not)

	Input infos for running with default values:
		- \-o <path/to/outputDir>, --out <path/to/outputDir>
					Name of output directory
		- \-c, --compress
						gzip output file

"""


##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import dict2txt, loadInList, dictDict2txt, printCol, relativeToAbsolutePath, existant_file, directory

## Python modules
import argparse
from time import localtime, strftime
import gzip

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='01/03/2017'

##################################################
## Functions

##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''Programme extract U from Tab file and build 3 tab:
	- prefilter: without N for all samples
	- withoutN: not missing data allow for all samples
	- withoutNandR: only line with one or more SNP''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-t', '--tab', metavar="<path/to/tabFileDir>",type = directory, required=True, dest = 'filesDir', help = 'path to file tab (gzip or not) ')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-o', '--out', metavar="<path/to/outputDir>",type = directory, default="./", dest = 'pathOut', help = 'Name of output directory (must exist)')
	files.add_argument('-c', '--compress',action ='store_true', dest = 'compress', help = 'gzip output file')


	# Check parameters
	args = parser.parse_args()


	#Welcome message
	print("#################################################################")
	print("#              Welcome on "+__file__+" (Version " + version + ")               #")
	print("#################################################################\n\n")

	pathFiles = args.filesDir
	pathFilesOut = args.pathOut
	compress = args.compress

	# resume value to user
	print(" - Intput Info:")
	print("\t - TAB files are in : %s" % pathFiles.pathDirectory)
	print("\t - %s TAB file count" % len(pathFiles.lsExtInDirToList(["tab","gz"])))

	print(" - Output Info:")
	if compress:
		print("\t - Output tab will be gzip")
	print("\t - Output tab created into directory:  %s\n\n" % pathFilesOut.pathDirectory)



	for tabFile in pathFiles.lsExtInDirToList(["tab","gz"]):

		# récupère le nom du fichier tab In
		basename = tabFile.split("/")[-1].split(".")[0]
		print(basename)
		outFileNamePrefilterN = pathFilesOut.pathDirectory+basename+"_prefilterN.tab"
		outFileNameWithoutN = pathFilesOut.pathDirectory+basename+"_withoutN.tab"
		outFileNamewithoutNandR = pathFilesOut.pathDirectory+basename+"_withoutNandR.tab"

		if compress:
			outFileNamePrefilterNFile = gzip.open(outFileNamePrefilterN+".gz", "wb")
			outFileNameWithoutNFile = gzip.open(outFileNameWithoutN+".gz", "wb")
			outFileNamewithoutNandRFile = gzip.open(outFileNamewithoutNandR+".gz", "wb")
		else:
			outFileNamePrefilterNFile = open(outFileNamePrefilterN, "wb")
			outFileNameWithoutNFile = open(outFileNameWithoutN, "wb")
			outFileNamewithoutNandRFile = open(outFileNamewithoutNandR, "wb")

		if '.gz' in tabFile:
			tabFileIn = gzip.open(tabFile, "rb")
		else:
			tabFileIn = open(tabFile, "rb")

		header = tabFileIn.readline().decode("utf-8")
		print(header)
		samples = header.rstrip().split("\t")[3:]
		nbSample = len(samples)


			outFileNamePrefilterNFile.write(bytes(header,"utf-8"))
			outFileNameWithoutNFile.write(bytes(header,"utf-8"))
			outFileNamewithoutNandRFile.write(bytes(header,"utf-8"))

		nbtotal, withoutN, prefilter,withoutNandR = 0, 0, 0, 0
		for line in tabFileIn:
			line = line.decode("utf-8")
			chrom, pos, ref = line.rstrip().split("\t")[:3]
			genotypes = line.rstrip().split("\t")[3:]

			nbN = genotypes.count("N")
			nbRef = genotypes.count(ref)
			if nbN != nbSample:
				if nbN == 0:
					withoutN+=1
					outFileNameWithoutNFile.write(bytes(line,"utf-8"))

					if nbRef!= nbSample:
						withoutNandR+=1
						outFileNamewithoutNandRFile.write(bytes(line,"utf-8"))
				else:
					prefilter+=1
					outFileNamePrefilterNFile.write(bytes(line,"utf-8"))
			nbtotal+=1
		outFileNamePrefilterNFile.close()
		outFileNameWithoutNFile.close()
		outFileNamewithoutNandRFile.close()
		tabFileIn.close()

		print("NBligne total: "+str(nbtotal))
		print("NBligne prefilter: "+str(prefilter))
		print("NBligne withoutN: "+str(withoutN))
		print("NBligne withoutNandR: "+str(withoutNandR))



	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
