#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
## @package splitFasta.py.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath, existant_file, fasta2dict

## Python modules
import argparse
from time import localtime, strftime

## BIO Python modules
from Bio import SeqIO


##################################################
## Variables Globales
version="0.1"
VERSION_DATE='16/06/2016'
debug="False"
#debug="True"


##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='splitFasta.py', description='''This Programme extract Fasta Seq with liste keep''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display splitFasta.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	files = parser.add_argument_group('Input info for running')
	files.add_argument('-f', '--fasta', metavar="<filename>", type=existant_file, required=True, dest = 'fastaFile', help = 'fasta files to split')
	files.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'pathFileOut', help = 'Name of output file directory')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#        Welcome in splitFasta (Version " + version + ")          #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	fastaFile = relativeToAbsolutePath(args.fastaFile)
	pathFileOut = relativeToAbsolutePath(args.pathFileOut)

	print("\t - Input Path is: %s" % pathFileOut.pathDirectory)
	print("\t - fasta file is : %s" % fastaFile)


	dicoFasta = fasta2dict(fastaFile)

	for name, sequence in dicoFasta.items():
		with open(pathFileOut.pathDirectory+name+".fasta", "w") as output_handle:
			SeqIO.write(sequence,output_handle, "fasta")

	#print("\n\nExecution summary:")

	#print("  - Outputting \n\
	#Il y a au final %i Sequences gardées sur les %i initial\n\
	#les sequences sont ajoutées dans le fichier %s" %(nbKeep,nbTotal,outputfilename))
	#print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
