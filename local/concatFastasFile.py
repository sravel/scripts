#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
# @package concatFastasFile.py
# @author Sebastien Ravel

"""
	The concatFastasFile script
	===========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 11/07/2016
	:version: 0.1

	Script description
	------------------

	This Programme concat multiple fasta files with same sequences name's into uniq file

	Example
	-------

	>>> concatFastasFile.py -d NT_ALIGN/ -o 2241Ortho-82souches.fasta

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display concatFastasFile.py version number and exit

	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path to directory fasta files (fasta","fas","fa","fna")
		- \-o <filename>, --out <filename>
						Name of output file

"""



##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath, directory, concatFastasFiles

# Python modules
import argparse
from time import localtime, strftime

## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='11/07/2016'
debug="False"
#debug="True"

##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='concatFastasFile.py', description='''This Programme concat multiple fasta files with same sequences name's into uniq file''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display concatFastasFile.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-d', '--directory', metavar="<path/to/directory>",type=directory, required=True, dest = 'fastaFileDir', help = 'path to directory fasta files ("fasta","fas","fa","fna")')
	filesReq.add_argument('-o', '--out', metavar="<filename>", required=True, dest = 'paramoutfile', help = 'Name of output file')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#        Welcome in concatFastasFile (Version " + version + ")          #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les arguments
	pathFastaFile = args.fastaFileDir
	outputfilename = relativeToAbsolutePath(args.paramoutfile)

	# resume value to user
	print(" - Intput Info:")
	print("\t - Fasta were in directory: %s" % pathFastaFile.pathDirectory)

	print(" - Output Info:")
	print("\t - Output file fasta is: %s" % outputfilename)


	nbfile = len(pathFastaFile.lsExtInDirToList(["fasta","fas","fa","fna"]))

	dico_concate = concatFastasFiles(pathFastaFile.pathDirectory)
	output_handle = open(outputfilename, "w")
	for ID, sequence in dico_concate.items():

		record = SeqRecord(sequence,id=ID,name=ID, description="")
		SeqIO.write(record,output_handle, "fasta")





	print("\n\nExecution summary:")

	print("  - Outputting \n\
	Il y a au final %i fichiers concaténer\n\
	les sequences sont ajouter dans le fichier %s\n" %(nbfile,outputfilename))
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
