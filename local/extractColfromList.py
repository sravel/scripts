#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
# @package extractColfromList.py
# @author Sebastien Ravel Lea Picard

"""
	The extractColfromList script
	=============================
	:author: Sebastien Ravel, Lea Picard
	:contact: sebastien.ravel@cirad.fr
	:date: 08/07/2016
	:version: 0.1

	Script description
	------------------

	This Program takes a list of IDs and extract the columns with corresponding IDs from a table

	Example
	-------

	>>> extractColfromList.py -pi NT_ALIGN/

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display extractColfromList.py version number and exit

	Input mandatory infos for running:
		- \-ti <filename>, --tablein <filename>
						Name of table file in
		- \-l <filename>, --listID <filename>
						File with IDs to be kept

	Input infos for running with default values:
		- \-to <filename>, --tableout <filename>
						Name of table file out (default = "tablein"_extractedIDs.tab)

"""


##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import loadInList, relativeToAbsolutePath, existant_file

from time import localtime, strftime
## Python modules
import argparse
import gzip

##################################################
## Variables Globales
version="0.2"
VERSION_DATE='20/03/2017'
#debug="False"
#debug="True"

##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='extractColfromList.py', description='''This Program takes a list of IDs and extract the columns with corresponding IDs from a table''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display extractColfromList version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-ti', '--tablein', metavar="<filename>", type=existant_file, required=True, dest = 'tableFile', help = 'Name of table file in')
	filesreq.add_argument('-l', '--listID', metavar="<filename>", type=existant_file, required=True, dest = 'IDlist', help = 'File with IDs to be kept')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-to', '--tableout', metavar="<filename>",default="", required=False, dest = 'tableFileOut', help = 'Name of table file out (default tablein_extractedIDs.tab)')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#            Welcome in extractColfromList (Version " + version + ")              #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	#get arguments
	tableFile = relativeToAbsolutePath(args.tableFile)
	tableFileOut = args.tableFileOut
	IDlist = args.IDlist

	if tableFileOut == "":
		tableFileOut=tableFile.split(".")[0]+"_extractedIDs.tab"

	#loading IDs to be kept in a list
	listNameKeep = loadInList(IDlist)

	if ".gz" in tableFile:
		fichier = gzip.open(tableFile,"rb")
	else:
		fichier = open(tableFile,"rb")


	#loading column IDs in a list
	header = fichier.readline().decode("utf-8").rstrip().split("\t")

	indice = 0
	listIndiceKeep = []
	for colName in header:
		if colName in listNameKeep:
			#print(colName, indice)
			#get position for each ID to be kept, +1 for cut command (sh col1 = python col0)
			listIndiceKeep.append(str(indice+1))
		indice += 1

	txtListIndiceKeep = ",".join(listIndiceKeep)
	#print(txtListIndiceKeep)

	print("\n\ncut -f"+txtListIndiceKeep+" "+tableFile+" > "+tableFileOut)
	os.system("cut -f"+txtListIndiceKeep+" "+tableFile+" > "+tableFileOut)

	print("\n\nExecution summary:")

	print("  - Outputting \n\
	Le ficher est créé dans %s \n" %(tableFileOut))

	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
