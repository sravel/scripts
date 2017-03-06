#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package mergeTabFiles.py
# @author Sebastien Ravel

"""
	The mergeTabFiles script
	=========================

	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 02/03/2017
	:version: 0.1

	Script description
	------------------

	This Programme merged files (compress or not) into one file (compress or not)

	Example
	-------

	>>> mergeTabFiles.py -f /work/ali/89_Mo/ -o ./mergeTabFile.tab.gz -e tab gz -m N

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display compress.py version number and exit

	Input mandatory infos for running:
		- \-f <path/to/fileDir>, --files <path/to/fileDir>
							Path to files to be merged (compress with gzip or not)
		- \-o <path/to/outputFileName>, --out <path/to/outputFileName>
						Name of output file merged if .gz in name autocompress file
		- \-e <str> [<str> ...], --extension <str> [<str> ...]
						loop to extension files (multiple format allows, eg:
						-e tab gz for tab and gz files)

	Input infos for running with default values:
		- \-m <string>, --missing <string>
						remplace missing data with input value
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
import pandas
import gzip
from functools import reduce


##################################################
## Variables Globales
version="0.1"
VERSION_DATE='27/01/2017'


##################################################
## Functions

def make_df(filename):
	if ".gz" in filename:
		print(filename)
		return pandas.read_csv(filename, sep="\t", compression="gzip")
	else:
		print(filename)
		return pandas.read_csv(filename, sep="\t")

def join_dfs(ldf, rdf):
	return ldf.merge(rdf, how='outer', left_on=['Chromosome', 'Position','Reference'], right_on=['Chromosome', 'Position','Reference'])


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme merged files (compress or not) into one file (compress or not)''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--files', metavar="<path/to/fileDir>",type = directory, required=True, dest = 'filesDir', help = 'Path to files to be merged (compress with gzip or not)')
	filesreq.add_argument('-o', '--out', metavar="<path/to/outputFileName>", required=True, dest = 'outFile', help = 'Name of output file merged if .gz in name autocompress file')
	filesreq.add_argument('-e', '--extension', metavar="<str>",type = str, nargs='+', required=True, dest = 'extParam', help = 'loop to extension files (multiple format allows, eg: -e tab gz for tab and gz files)')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-m', '--missing', metavar="<string>",type=str, default="N", dest = 'missingParam', help = 'remplace missing data with input value')


	# Check parameters
	args = parser.parse_args()


	#Welcome message
	print("#################################################################")
	print("#              Welcome on mergeFiles (Version " + version + ")               #")
	print("#################################################################\n\n")

	pathFiles = args.filesDir
	extParam = args.extParam
	missingParam = args.missingParam

	outFile = args.outFile


	# resume value to user
	print(" - Intput Info:")
	print("\t - Files are in : %s" % pathFiles.pathDirectory)
	print("\t - Files extension merged : '%s'" % ",".join(extParam))
	print("\t - %s files count" % len(pathFiles.lsExtInDirToList(extParam)))
	print("\t - Missing data remplace by: %s" % missingParam)

	print(" - Output Info:")
	print("\t - Output Merge tab created:  %s\n\n" % outFile)


	fileList = []
	for fileIn in pathFiles.lsExtInDirToList(extParam):
		fileList.append(fileIn)

	#print("Files are\n\t%s \n\n" % "\n\t".join(fileList))

	dfs = [make_df(filename) for filename in fileList]

	final_df = reduce(join_dfs, dfs)
	print("Head of final merged: (NaN remplace by %s in outFile)\n" % missingParam)
	print(final_df.head())

	# sauvegarder la dataframe dans un fichier
	if ".gz" in outFile:
		final_df.to_csv(outFile, sep="\t", index=False, compression="gzip", na_rep=missingParam)
	else:
		final_df.to_csv(outFile, sep="\t", index=False, na_rep=missingParam)




	print("\n\n#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
