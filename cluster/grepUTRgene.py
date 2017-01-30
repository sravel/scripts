#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package grepUTRgene.py
# @author Sebastien Ravel

"""
	The grepUTRgene script
	=========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 08/11/2016
	:version: 0.1

	Script description
	------------------

	This Programme build fasta of listed genes with UTR size

	Example
	-------

	>>> grepUTRgene.py -c chromosome.fasta -r listGenesID.txt -g gffFile.gff -o genesUTR.fasta -i ID [-u 3000]

	Help Programm
	-------------

	optional arguments:
		- \-h, --help            show this help message and exit
		- \-v, --version         display ./grepUTRgene.py version number and exit
		- \-dd, --debug          enter verbose/debug mode

	Input mandatory infos for running:
		- \-c <path/to/chromosomeFile>, --chro <path/to/chromosomeFile>
						path to chromosomeFile files
		- \-l <path/to/listGene>, --listgene <path/to/listGene>
						path to list of gene to extract
		- \-g <path/to/gff>, --gff <path/to/gff>
						path to gff file
		- \-o <path/to/outputName>, --out <path/to/outputName>
						Name of output file
		- \-i <string>, --id <string>
						String of you id tag in gff (like "proteinId",
						ID",...)

	Input infos for running with default values:
		- \-u <int>, --utr <int>
						len of UTR keep (default = 3000)

"""

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath, directory, printCol, existant_file, loadInList, parseGFF

## Python modules
import argparse
from time import localtime, strftime

## BIO Python modules
## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet

from pyfaidx import Fasta

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='10/11/2016'

##################################################
## Functions

###################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme build fasta of listed genes with UTR size''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug',action='store_true', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-c', '--chro', metavar="<path/to/chromosomeFile>", type=existant_file, required=True, dest = 'chromosomeFile', help = 'path to chromosomeFile files')
	filesReq.add_argument('-l', '--listgene', metavar="<path/to/listGene>", type=existant_file, required=True, dest = 'listGeneFile', help = 'path to list of gene to extract')
	filesReq.add_argument('-g', '--gff', metavar="<path/to/gff>", type=existant_file, required=True, dest = 'gffFile', help = 'path to gff file')
	filesReq.add_argument('-o', '--out', metavar="<path/to/outputName>", required=True, dest = 'outFile', help = 'Name of output file')
	filesReq.add_argument('-i', '--id', metavar="<string>", required=True, dest = 'idTag', help = 'String of you id tag in gff (like "proteinId", "ID",...) ')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-u', '--utr', metavar="<int>", type = int, default=3000,required=False, dest = 'UTRlen', help = 'len of UTR keep (default = 3000)')
	#files.add_argument('-th', '--thread', metavar="<int>",type = int, default=4, required=False, dest = 'nbThreads', help = 'number of threads for mapping (default = 4)')


	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#           Welcome in %s (Version %s)            #" %(__file__, version))
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les arguments
	chromosomeFile = args.chromosomeFile
	listGeneFile = args.listGeneFile
	gffFile = args.gffFile
	outFile = args.outFile
	UTRlen = args.UTRlen


	# resume value to user
	print(" - Intput Info:")
	print("\t - Chromosome file is: %s" % chromosomeFile)
	print("\t - List of gene file is: %s" % listGeneFile)
	print("\t - GFF file is : %s" % gffFile)
	print("\t - UTRlen is : %s" % UTRlen)
	print("\t - TAG is : %s" % args.idTag)

	print(" - Output Info:")
	print("\t - Output file with gene+UTR is:  %s" % outFile)


	# chargement du fasta des chromosome en mémoire
	sequencesChrom = Fasta(chromosomeFile)
	#print(sequencesChrom.keys())

	# chargement de la list des gènes
	keepGeneList = loadInList(listGeneFile)
	print("\t - There is %s genes" % len(keepGeneList))

	# ouverture du fichier de sortie:
	with open(outFile, "w") as outputFile:
		# parse GFF pour avoir les positions
		objGFF = parseGFF(gffFile)
		for record in objGFF.parseGFF3():
			chromosome = record.seqid
			if record.type == "gene":
				try:
					geneName = record.attributes[args.idTag]
				#print(geneName)
				except Exception as e:
					pass

				if geneName in keepGeneList:
					start, stop = record.start, record.end
					UTRstart, UTRstop = int(record.start)-UTRlen, int(record.end)+UTRlen
					#print("%s\t%s\t%s" %(geneName,UTRstart, UTRstop))

					if UTRstart <= 0:
						UTRstart = 1

					if start != 0:
						seq = sequencesChrom[chromosome][int(UTRstart)-1:int(UTRstop)].seq
					elif start == 0:
						seq = sequencesChrom[chromosome][int(UTRstart):int(UTRstop)].seq
					lenSeq = len(seq)
					des = "%s:%s..%s UTR-%s lenght=%s" % (chromosome,UTRstart,UTRstop, UTRlen, lenSeq)

					record = SeqRecord(Seq(seq),id=geneName,name=geneName, description=des)
					SeqIO.write(record,outputFile, "fasta")



	#print("\n - Execution summary:")

	#print("\n  You want run MutilmappingX for %s fasta,\
 #The script are created all fasta-MutilmappingX.sh for all fasta into %s,\n\
 #For run all sub-script in qsub, %s was created.\n\
 #It lunch programm with job array and run %s job max:\n" %(count-1,outputSHDir,SGENameFile, args.nbJobValue))
	#printCol.green("\tqsub %s" % SGENameFile)
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
