#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
## @package extractSeqFasta.py.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys
sys.path.insert(1,'../modules/')
from MODULES_SEB import extractInverseListFromFasta, extractListFromFasta, relativeToAbsolutePath

## Python modules
import argparse
from time import localtime, strftime

## BIO Python modules
from Bio import SeqIO


##################################################
## Variables Globales
version="0.1"
VERSION_DATE='04/03/2015'
debug="False"
#debug="True"


##################################################
## Functions
def checkParameters (arg_list):
	# Check input related options
	if (not arg_list.fastaFile):
		print ('Error: No input file defined via option -f/--fasta !' + "\n")
		parser.print_help()
		exit()

##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='extractSeqFasta.py', description='''This Programme extract Fasta Seq with liste keep''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display extractSeqFasta.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	files = parser.add_argument_group('Input info for running')
	files.add_argument('-f', '--fasta', metavar="<filename>", dest = 'fastaFile', help = 'fasta files')
	files.add_argument('-l', '--list', metavar="<filename>", dest = 'listFile', help = 'list files with keep name or not keep')
	files.add_argument('-o', '--out', metavar="<filename>", dest = 'paramoutfile', help = 'Name of output file')
	files.add_argument('-k', '--keep', metavar="<yes/y/no/n>", dest = 'keepValue',choices = ["yes","y","no","n"], help = 'choise keep (y/yes) or not keep (n/no) sequences in list file')

	# Check parameters
	args = parser.parse_args()
	checkParameters(args)

	#Welcome message
	print("#################################################################")
	print("#        Welcome in extractSeqFasta (Version " + version + ")          #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	fastaFile = relativeToAbsolutePath(args.fastaFile)
	outputfilename = relativeToAbsolutePath(args.paramoutfile)
	listFile = relativeToAbsolutePath(args.listFile)
	keepValue = args.keepValue

	output_handle = open(outputfilename, "w")


	if keepValue in ["no","n"]:
		dico_keep, nbTotal = extractInverseListFromFasta(fastaFile, listFile)
	elif keepValue in ["yes","y"]:
		dico_keep, nbTotal = extractListFromFasta(fastaFile, listFile)

	nbKeep = len(dico_keep.keys())
	for geneId, sequence in dico_keep.items():
		SeqIO.write(sequence,output_handle, "fasta")

	output_handle.close()
	print("\n\nExecution summary:")

	print("  - Outputting \n\
	Il y a au final %i Sequences garder sur les %i initial\n\
	les sequences sont ajouter dans le fichier %s" %(nbKeep,nbTotal,outputfilename))
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
