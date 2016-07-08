#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package buildSNPtoFasta.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath, existant_file, fasta2dict, directory, parseGFF

## Python modules
import argparse
from time import localtime, strftime

## BIO Python modules
from Bio import SeqIO


##################################################
## Variables Globales
version="0.1"
VERSION_DATE='08/07/2016'
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
	parser = argparse.ArgumentParser(prog='buildSNPtoFasta.py', description='''This Programme parse GFF, TAB to add SNP in sequences''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display buildSNPtoFasta.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	files = parser.add_argument_group('Input info for running')
	files.add_argument('-g', '--gff', metavar="<filename>", type=existant_file, required=True, dest = 'gffFile', help = 'gff files with annotation')
	#files.add_argument('-f', '--fasta', metavar="<filename>", type=existant_file, required=True, dest = 'fastaFile', help = 'fasta files to split')
	#files.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'pathOut', help = 'Name of output file directory')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#        Welcome in buildSNPtoFasta (Version " + version + ")          #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	gffFile = relativeToAbsolutePath(args.gffFile)
	#fastaFile = relativeToAbsolutePath(args.fastaFile)
	#pathFileOut = args.pathOut

	print("\t - Input GFF is: %s" % gffFile)
	#print("\t - Input Path is: %s" % pathFileOut.pathDirectory)
	#print("\t - fasta file is : %s" % fastaFile)


	objGFF = parseGFF(gffFile)
	recordCount = 0

	for record in objGFF.parseGFF3():
		#Apply filter, if any
		#if args.filter_type and record.type != args.filter_type:
			#continue
		#Print record if specified by the user
		if record.type == "CDS":
			print(record)
		#Access attributes like this: my_strand = record.strand
		recordCount += 1
	print("Total records: %d" % recordCount)
