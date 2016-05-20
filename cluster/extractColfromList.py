#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
## @package extractColfromList.py
# @author Lea Picard

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import loadInList, relativeToAbsolutePath

from time import localtime, strftime
## Python modules
import argparse

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='04/01/2016'
#debug="False"
#debug="True"

##################################################
## Main code
##################################################
if __name__ == "__main__":


	# Parameters recovery
	parser = argparse.ArgumentParser(prog='extractColfromList.py', description='''This Program takes a list of IDs and extract the columns with corresponding IDs from a table''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display extractColfromList version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	files = parser.add_argument_group('Input info for running')
	files.add_argument('-ti', '--tablein', metavar="<filename>", dest = 'tableFile', help = 'Name of table file in')
	files.add_argument('-l', '--listID', metavar="<filename>", dest = 'IDlist', help = 'File with IDs to be kept')
	files.add_argument('-to', '--tableout', metavar="<filename>", dest = 'tableFileOut', help = 'Name of table file out (default tablein_extractedIDs.tab)')

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

	fichier = open(tableFile,"r")

	#loading column IDs in a list
	header = fichier.readline().rstrip().split("\t")

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

	os.system("cut -f"+txtListIndiceKeep+" "+tableFile+" > "+tableFileOut)

	print("\n\nExecution summary:")

	print("  - Outputting \n\
	Le ficher est créé dans %s \n" %(tableFileOut))

	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
