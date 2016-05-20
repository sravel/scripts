#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
## @package extractSeqFastaFromLen.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
## Python modules
import argparse
from time import localtime, strftime
## BIO Python modules
from Bio import SeqIO

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='31/03/2015'
debug="False"
#debug="True"


##################################################
## Functions
def fasta2dict(filename):
	"""
	Function that take a file name (fasta), and return a dictionnary of sequence

	:param filename: a fasta file
	:type filename: file in fasta format
	:rtype: record_dict()
	:return: dict() - dictionnary with keys are Id and value SeqRecord() fields
	:requires: this function require ## BIO Python modules: (from Bio import SeqIO,\\n
															 from Bio.SeqRecord import SeqRecord \\n
															 from Bio.Seq import Seq \\n
															 from Bio.Alphabet import SingleLetterAlphabet)

	Example:
		>>> filename = "sequence.fasta"
		>>> fasta2dict(filename)
		{">Seq1":"SeqRecord()"}
	"""

	# chargement du fasta des MGG en mémoire
	handle = open(filename, "rU")
	record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
	handle.close()
	return record_dict

def lenSeq2dict(filename):
	"""
	Function that take a file name (fasta), and return a dictionnary with length of sequence

	:param filename: a fasta file
	:type filename: file in fasta format
	:rtype: record_dict()
	:return: dict() - contain length of sequences
	:requires: this function require fasta2dict(filename)

	Example:
		>>> filename = "sequence.fasta"
		>>> lenSeq2dict(filename)
		{">Seq1":20154}
	"""

	dicoLenMGG = {}
	record_dict = fasta2dict(filename)
	for gene in sorted(record_dict.keys()):
		if record_dict[gene].id not in dicoLenMGG:
			lenseq = len(record_dict[gene].seq)
			dicoLenMGG[gene]=int(lenseq)
	return dicoLenMGG

def relativeToAbsolutePath(relative):
	"""	Return the absolutPath

	:param relative: a string path
	:type relative: string
	:rtype: string()
	:return: absolutePath
	:warn: need subprocess::check_output

	Example:
		>>> print(relative)
			../test
		>>> pathDirectory = relativeToAbsolutePath(relative)
		>>> print(pathDirectory)
			/home/sebastien/test

	"""
	from subprocess import check_output
	if relative[0] != "/":			# The relative path is a relative path, ie do not starts with /
		command = "readlink -m "+relative
		absolutePath = subprocess.check_output(command, shell=True).decode("utf-8").rstrip()
		return absolutePath
	else:						# Relative is in fact an absolute path, send a warning
		absolutePath = relative;
		return absolutePath


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='extractSeqFastaFromLen.py', description='''This Programme extract Fasta Seq with liste keep''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display extractSeqFastaFromLen.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	files = parser.add_argument_group('Input info for running')
	files.add_argument('-f', '--fasta', metavar="<filename>", required=True, dest = 'fastaFile', help = 'fasta files')
	files.add_argument('-l', '--len', metavar="<int>", required=True, type = int, dest = 'lenSize', help = 'lensize cutoff')
	files.add_argument('-o', '--out', metavar="<filename>", required=True, dest = 'paramoutfile', help = 'Name of output file')
	files.add_argument('-k', '--keep', metavar="<g/greater/l/lower>", required=True, dest = 'keepValue',choices = ["g","greater","l","lower"], help = 'choice keep sequences size greater than -l (g/greater) or keep lower (l/lower)')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#        Welcome in extractSeqFastaFromLen (Version " + version + ")          #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	fastaFile = relativeToAbsolutePath(args.fastaFile)
	outputfilename = relativeToAbsolutePath(args.paramoutfile)
	lenSize = args.lenSize
	keepValue = args.keepValue
	output_handle = open(relativeToAbsolutePath(outputfilename), "w")

	dicoSize = lenSeq2dict(fastaFile)
	dicoFasta = fasta2dict(fastaFile)

	nbKeep=0
	nbTotal = len(dicoFasta.keys())

	for ID, lenSeq in dicoSize.items():
		if keepValue in ["g","greater"]:
			if lenSeq >= lenSize:
				sequence = dicoFasta[ID]
				SeqIO.write(sequence,output_handle, "fasta")
				nbKeep += 1

		elif keepValue in ["l","lower"]:
			if lenSeq <= lenSize:
				sequence = dicoFasta[ID]
				SeqIO.write(sequence,output_handle, "fasta")
				nbKeep += 1



	print("\n\nExecution summary:")

	print("  - Outputting \n\
	Il y a au final %i Sequences garder sur les %i initial\n\
	les sequences sont ajouter dans le fichier %s" %(nbKeep,nbTotal,outputfilename))
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
