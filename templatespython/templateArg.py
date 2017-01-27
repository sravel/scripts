#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package scriptName.py
# @author Sebastien Ravel

"""
	The scriptName script
	===========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 08/07/2016
	:version: 0.1

	Script description
	------------------

	This Programme ....

	Example
	-------

	>>> scriptName.py -d asie_2015_480mlg -c CLUMPAK/ -l asie_480mlg_STRUC.txt

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display scriptName.py version number and exit

	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of result structure
		- \-c <path/to/directory/clumpak>, --clumpak <path/to/directory/clumpak>
						path of clumpak directory
		- \-l <filename>, --label <filename>
						File with LABEL, first column name, second top label
						info

	Input infos for running with default values:
		- \-dp <filename>, --drawparams <filename>
						Check your own drawparams file
		- \-co <filename>, --color <filename>
						File with colors (default 15 color max)

"""


##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import directory, relativeToAbsolutePath, existant_file, printCol

## Python modules
import argparse
from time import localtime, strftime

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='27/01/2017'


##################################################
## Functions


##################################################
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

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory', metavar="<path/to/directory>",type = directory, required=True, dest = 'dirPath', help = 'path of result structure')
	filesreq.add_argument('-c', '--clumpak', metavar="<path/to/directory/clumpak>",type = directory, required=True, dest = 'dirPathClumpak', help = 'path of clumpak directory')
	filesreq.add_argument('-l', '--label', metavar="<filename>",type=existant_file, required=True, dest = 'labelFileParam', help = 'File with LABEL, first column name, second top label info')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-dp', '--drawparams', metavar="<filename>",type=existant_file, required=False, dest = 'drawparamsParam', help = 'Check your own drawparams file')
	files.add_argument('-co', '--color', metavar="<filename>",type=existant_file, required=False, dest = 'colorParam', help = 'File with colors (default 15 color max)')


	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#           Welcome in %s (Version %s)            #" %(__file__, version))
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les infos passer en argument
	workingObjDir = args.dirPath
	clumpakObjDir = args.dirPathClumpak
	labelFileParam = relativeToAbsolutePath(args.labelFileParam)


	print("\n\nExecution summary:")
	print("  - Outputting \n\
\t- Files created in %s." % nameDirectoryOutput)

	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
