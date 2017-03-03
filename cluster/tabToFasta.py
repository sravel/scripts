#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package tabToFasta.py
# @author Sebastien Ravel

"""
	The tabToFasta script
	===========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 02/03/2017
	:version: 0.1

	Script description
	------------------

	This Programme build fasta from tab file Adapt the RAM according to the size of the array (eg for 33127958 line and 48 individuals uses 34GB of RAM)

	Example
	-------

	>>> tabToFasta.py t 48Mo_table_withoutN.tab

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display tabToFasta.py version number and exit

	Input mandatory infos for running:
		- \-t <filename>, --tab <filename>
						Tab file

	Input infos for running with default values:
		- \-f <filename>, --fasta <filename>
						fasta Out (default = TabBasename.fasta)
		- \-l <filename>, --list <filename>
						Change Individual ID with custom ID provied table (tab
						with fisrt col ID, second col custom ID)
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
from MODULES_SEB import directory, relativeToAbsolutePath, existant_file, printCol, dict2txt, loadInDictCol

## Python modules
import argparse
from time import localtime, strftime
import gzip, shutil

## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='27/01/2017'


##################################################
## Functions
def transpose_matrix(matrix):
	mt = list( map(list, zip(*matrix)))
	return mt

		#zip(*matrix)
			#return matrix

##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme build fasta from tab file Adapt the RAM according to the size of the array (eg for 33127958 line and 48 individuals uses 34GB of RAM)''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug', action='store_true', dest='debug', help='enter verbose/debug mode')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-t', '--tab', metavar="<filename>",type=existant_file, required=True, dest = 'tabFileParam', help = 'Tab file (gzip or not)')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-f', '--fasta', metavar="<filename>", required=False, dest = 'outFastaParam', help = 'fasta Out (default = TabBasename.fasta)')
	files.add_argument('-l', '--list', metavar="<filename>", type=existant_file, required=False, dest = 'IDParam', help = 'Change Individual ID with custom ID provied table (tab with fisrt col ID, second col custom ID), include Reference ID')
	files.add_argument('-c', '--compress',action ='store_true', dest = 'compress', help = 'gzip output file ')



	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#           Welcome in %s (Version %s)            #" %(__file__, version))
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les infos passer en argument
	tabFileParam = args.tabFileParam
	outFastaParam = args.outFastaParam
	IDParam = args.IDParam
	compress = args.compress

	basename = tabFileParam.split("/")[-1].split(".")[0]
	print(basename)

	if outFastaParam == None:
		outFastaParam= relativeToAbsolutePath(basename+".fasta")

	# resume value to user
	print(" - Intput Info:")
	print("\t - TAB files is : %s" % tabFileParam)
	if IDParam != None:
		print("\t - Change Individual ID with custom ID provied table : %s" % IDParam)
		dicoCustomID = loadInDictCol(IDParam,0,1)

	print(" - Output Info:")
	if compress:
		print("\t - Output fasta will be gzip")
	print("\t - Output fasta is:  %s\n\n" % outFastaParam)


	if ".gz" in tabFileParam:
		tabFileIn =  gzip.open(tabFileParam, "rb")
	else:
		tabFileIn =  open(tabFileParam, "rb")

	outFile = open(outFastaParam, "w")

	matrice = [line.decode("utf-8").rstrip().split("\t") for line in tabFileIn]
	matriceT = transpose_matrix(matrice)
	matrice=""

	for listSample in matriceT[2:]:
		sample = listSample[0]
		seq = "".join(listSample[1:])

		if IDParam != None:
			sample = dicoCustomID[sample]
		record = SeqRecord(Seq(seq),id=sample,name=sample, description="")
		SeqIO.write(record,outFile, "fasta")


	tabFileIn.close()
	outFile.close()

	if compress:
		with open(outFastaParam, 'rb') as f_in:
			with gzip.open(outFastaParam+'.gz', 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)

	############### Other Methode

		#with open(tabFileParam, "r") as tabFileIn:

		#header = tabFileIn.readline()
		#samples = header.rstrip().split("\t")[2:]
		#nbSample = len(samples)
		#dicoSeq = {}
		#for sample in samples:
			#dicoSeq[sample] = ""

		#printCol.green('number of samples: %s\n' % nbSample)
		#printCol.green('samples: %s' % "\t".join(samples))

		#for line in tabFileIn:
			#genotypes = line.rstrip().split("\t")[2:]

			#i=0
			#for genotype in genotypes:
				#sample = samples[i]
				#dicoSeq[sample]+=genotype
				#i+=1
	##print(dict2txt(dicoSeq))
	#with open(outFastaParam, "w") as outFile:
		#for sample, seq in dicoSeq.items():
			#record = SeqRecord(Seq(seq),id=sample,name=sample, description="")
			#SeqIO.write(record,outFile, "fasta")


	print("\n\nExecution summary:")


	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
