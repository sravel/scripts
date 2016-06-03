#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
## @package countCDSStartStop.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import directory, sort_human, fasta2dict

## Python modules
import argparse, re
from time import localtime, strftime

## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='03/05/2016'
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
# 	start=time.clock()
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='countCDSStartStop.py', description='''This Programme count CDS with Start AND/OR Stop or NONE Both, in fasta files''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display countCDSStartStop version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-pi', '--pathin', metavar="<path/to/filein>", type = directory, required=True, dest = 'pathFileIn', help = 'path to fasta files in')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-f', '--fasta', metavar="<yes/y/no/n>", default="no", dest = 'fastaValue',choices = ["yes","y","no","n"], help = 'choice make (y/yes) or not (n/no) real fasta file (default=no)')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#        Welcome in countCDSStartStop (Version " + version + ")           #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	workingObjDir = args.pathFileIn

	print("Workink Directory: %s" % workingObjDir.pathDirectory)

	# ajoute à la variable current_dir le chemin ou est executer le script
	current_dir = os.path.dirname(os.path.abspath(__file__))+"/"

	startATG=0
	stop=0
	nbstart=0
	nbstop=0
	nbstartandstop=0
	nbEmpty=0
	other=0
	listCDSfiles = workingObjDir.lsExtInDirToList("codingseq")
	statFile = open(workingObjDir.pathDirectory+"statsCDSInfo.txt","w")
	listFile = open(workingObjDir.pathDirectory+"CDSfilelist.txt","w")

	statFile.write("fileName\tnbstartonly\tnbstoponly\tnbstartandstop\tnbEmpty\ttotal\tOther\n")

	for fileCDS in sorted(listCDSfiles):
		fileName = fileCDS.split("/")[-1].split(".")[0]
		listFile.write(fileCDS+"\n")

		if args.fastaValue in ["yes","y"]:
			output_file = open(workingObjDir.pathDirectory+fileName+".fasta", "w")

		record_dict = fasta2dict(fileCDS)
		for name in sorted(record_dict.keys(), key=sort_human):
			record = record_dict[name]
			#oldNumID = record.id
			#new_record_name = fileName+"_"+oldNumID
			#record.id = new_record_name
			#record.name = ""
			seq = record.seq.upper()
			firstCodon = seq[:3]
			endCodon = seq[-3:]

			# Test ATG start and Stop codons
			if str(firstCodon.upper()) in "ATG":
				startATG=1
			else:
				startATG=0
			if str(endCodon.upper()) in ["TAG","TAA","TGA"]:
				stop=1
			else:
				stop=0
			if startATG == 1 and stop == 1:
				nbstartandstop+=1
				record.description = 'start_stop'
			elif startATG == 1 and stop == 0:
				nbstart+=1
				record.description = 'start_-'
			elif startATG == 0 and stop == 1:
				nbstop+=1
				record.description = '-_stop'
			elif startATG == 0 and stop == 0:
				nbEmpty+=1
				record.description = '-_-'
			else:
				other+=1
			# write the whole thing out
			if args.fastaValue in ["yes","y"]:
				SeqIO.write(record.upper(), output_file, 'fasta')

		if nbstart !=0 or nbstartandstop !=0 or nbstop !=0:
			print(fileName)
			total = nbstart+nbstop+nbstartandstop+nbEmpty
			statFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(fileName,nbstart,nbstop, nbstartandstop,nbEmpty,total,other))
			nbstart=0
			nbstop=0
			nbstartandstop=0
			nbEmpty = 0
			other = 0
		if args.fastaValue in ["yes","y"]:
			output_file.close()

	statFile.close()
	listFile.close()



	print("\n\nExecution summary:")

	print("  - Outputting \n\
	- File with Stat résult output to : %s\n\
	- List of fasta files output to : %s\n" %(statFile.name,listFile.name))
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
