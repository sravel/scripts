#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package grepMotifFromAlignment.py
# @author Sebastien Ravel

"""
	The grepMotifFromAlignment script
	=================================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 08/07/2016
	:version: 0.1

	Script description
	------------------

	This Programme parse Aligment info to build motif table of SNP in gene's

	Example
	-------

	>>> grepMotifFromAlignment.py -d path/to/fasta -o filenameout

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display grepMotifFromAlignment.py version number and exit

	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path to directory fasta files
		- \-o <filename>, --out <filename>
						Name of output file
	Input infos for running with default values:


"""

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import directory, dictList2txt, dictDict2txt, dict2txt, relativeToAbsolutePath, existant_file

## Python modules
import argparse
from time import localtime, strftime
## BIO Python modules
from Bio import AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='04/03/2015'
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
	parser = argparse.ArgumentParser(prog='grepMotifFromAlignment.py', description='''This Programme parse Aligment info to build motif table of SNP in gene's''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display grepMotifFromAlignment.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory', metavar="<path/to/directory>",type=directory, required=True, dest = 'pathDirectory', help = 'path to directory fasta files')
	filesreq.add_argument('-o', '--out', metavar="<filename>", required=True, dest = 'paramoutfile', help = 'Name of output file')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-l', '--list', metavar="<filename>", default="ALL", dest = 'listKeepFile', help = 'File with Strain to keep (one per row)')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#        Welcome in grepMotifFromAlignment (Version " + version + ")          #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	pathDirectory = args.pathDirectory
	outputfilename = relativeToAbsolutePath(args.paramoutfile)

	if args.listKeepFile not in ["ALL"]:
		listKeepSouche = loadInList(existant_file(args.listKeepFile))
	else:
		listKeepSouche = ["ALL"]


	print("\t - Input pathDirectory is: %s" % pathDirectory)
	print("\t - Output file name is: %s" % outputfilename)
	print("\t - You want to keep strain:\n%s" % "\n".join(listKeepSouche))


	dicoOutputTxt = {}
	dicoSeqSNP = {}
	dicoFilenbSNP ={}
	dicoFileCountSNP ={}
	fileEmpty = 0
	listFileEmpty = []
	ctr = 1
	nbMotifTotal=0

	for filein in pathDirectory.listFiles:
		ctr += 1
		if ((ctr % 100 == 0) and (ctr != 0)) or (float(ctr) == len(pathDirectory.listFiles)):
			percent = (float(ctr)/float(len(pathDirectory.listFiles)))*100
			sys.stdout.write("\rProcessed up to %0.2f %%..." % percent)
			sys.stdout.flush()

		#print(filein)
		dicoSeqSNP = {}
		nbSNP = 0
		tableauSoucheName = []
		# lecture de l'alignement
		#alignment = AlignIO.read(open(resultataligment,"r"), "fasta")

		alignmentStart = AlignIO.read(open(filein,"r"), "fasta")
		# cree un nouvelle alignement avec que les souches voulus:
		keepListRecord = []

		for record in alignmentStart:
			if record.id not in listKeepSouche and "ALL" in listKeepSouche:
				listKeepSouche.append(record.id)
			#print(record.id)
			if record.id in listKeepSouche:
				keepListRecord.append(record)
				tableauSoucheName.append(record.id)
				if record.id not in dicoSeqSNP.keys():
					dicoSeqSNP[record.id] = ""
		alignment = MultipleSeqAlignment(keepListRecord)
		lenAlignement = int(alignment.get_alignment_length())
		#print(alignment)
		#print(tableauSoucheName)
		#print(len(tableauSoucheName))

		for indice in range(0,lenAlignement):
			tab = list(alignment[:,indice])
			#print(tab)
			nbO = tab.count(tab[0])
			nbA = tab.count("A")
			nbC = tab.count("C")
			nbT = tab.count("T")
			nbG = tab.count("G")
			nbN = tab.count("N")+tab.count("n")
			nbGap = tab.count("-")
			sommeACTG = nbA + nbC + nbT + nbG
			allcount = sommeACTG + nbN + nbGap
			if int(allcount) != len(alignment):						# test si total = nombre de souche
				print( sommeACTG, nbA , nbC , nbT, nbG,nbN, nbGap)
				print( tab)
				exit()
			if nbGap == 0 :
				if nbO != sommeACTG and nbN == 0:
					nbSNP+=1
					#print(indice)
					for lentabi in range(0,len(tab)):
						dicoSeqSNP[tableauSoucheName[lentabi]] += (tab[lentabi])

		nbSNPtotal=nbSNP
		if nbSNPtotal == 0:
			fileEmpty += 1
			listFileEmpty.append(filein)

		else:
			nbMotifTotal+=1
			listMotif = []
			for geneId, sequence in dicoSeqSNP.items():
				nbSNPtotal = (len(sequence))

				listMotif.append(sequence)

			nameMGG = filein.split("/")[-1].replace("_Orthologue_macse_NT.fasta","")
			if nameMGG not in dicoFileCountSNP.keys():
				dicoFileCountSNP[nameMGG] = {"NBSNP":nbSNPtotal,
											 "lenAlign":lenAlignement}


			if nbSNPtotal not in dicoFilenbSNP.keys():
				dicoFilenbSNP[nbSNPtotal] = 1
			else:
				dicoFilenbSNP[nbSNPtotal] +=1

			#print(nbSNPtotal)

			dicoCompteMotif = {k: listMotif.count(k) for k in set(listMotif)}
			#print(dict2txt(dicoCompteMotif))

			dicoTranslateMotif2Code = {}
			code = 10
			for motifUniq in dicoCompteMotif.keys():
				dicoTranslateMotif2Code[motifUniq] = code
				code+=1


			for geneId, sequence in dicoSeqSNP.items():
				codeSeq = dicoTranslateMotif2Code[sequence]
				if geneId not in dicoOutputTxt.keys():
					dicoOutputTxt[geneId] = [str(codeSeq)]
				else:
					dicoOutputTxt[geneId].append(str(codeSeq))

	output_handle = open(outputfilename, "w")
	#print(dictList2txt(dicoOutputTxt))
	outputTxt = ""
	#for key in sorted(dicoOutputTxt.keys()):
	for key in listKeepSouche:
		value = "\t".join(dicoOutputTxt[key])
		outputTxt += "%s\t%s\n" % (str(key),str(value))
	output_handle.write(outputTxt)

	outputListEmpty = open(basename+"_outputListEmpty.txt", "w")
	for fileEmptyName in listFileEmpty:
		outputListEmpty.write(fileEmptyName+"\n")


	with open(basename+"_LenAlign_nbSNP.txt","w") as output1:
		txt1 = dictDict2txt(dicoFileCountSNP)
		output1.write(txt1)


	with open(basename+"_nbSNPallFile.txt","w") as output2:
		txt1 = dict2txt(dicoFilenbSNP)
		output2.write(txt1)


	print("\n\nExecution summary:")

	print("  - Outputting \n\
	Il y a au final %i Motif dans tout les MGG\n\
	Il y a %i fichiers vides\n\
	les sequences sont ajouter dans le fichier %s\n\
	la liste des fichiers vides est dans le fichier outputListEmpty.txt" %(nbMotifTotal,fileEmpty,outputfilename))
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
