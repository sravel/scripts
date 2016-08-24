#!/usr/bin/python3.5
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
		listKeepFile = loadInList(existant_file(args.listKeepFile))
	else:
		listKeepFile = "ALL"


	print("\t - Input pathDirectory is: %s" % pathDirectory)
	print("\t - Output file name is: %s" % outputfilename)
	print("\t - You want to keep strain:\n%s" % "\n".join(listKeepFile))


	dicoOutputTxt = {}
	dicoSeqSNP = {}
	dicoFilenbSNP ={}
	dicoFileCountSNP ={}
	fileEmpty = 0
	listFileEmpty = []

	for filein in pathDirectory.listFiles:
		print(filein)
		dicoSeqSNP = {}
		nbSNP = 0
		tableauSoucheName = []
		# lecture de l'alignement
		#alignment = AlignIO.read(open(resultataligment,"r"), "fasta")
		#try:
		alignment = AlignIO.read(open(filein,"r"), "fasta")
		#print(alignment)
		lenAlignement = int(alignment.get_alignment_length())
		for record in alignment:
			tableauSoucheName.append(record.id)
			if record.id not in dicoSeqSNP.keys():
				dicoSeqSNP[record.id] = ""
			#else:
				#print("Exists")
		#print(tableauSoucheName)

		for indice in range(0,lenAlignement):
			tab = list(alignment[:,indice])
			nbO = tab.count(tab[0])
			nbA = tab.count("A")
			nbC = tab.count("C")
			nbT = tab.count("T")
			nbG = tab.count("G")
			nbN = tab.count("N")+tab.count("n")
			nbGap = tab.count("-")
			sommeACTG = nbA + nbC + nbT + nbG
			allcount = sommeACTG + nbN + nbGap
			if int(allcount) != len(tableauSoucheName):
				print( sommeACTG, nbA , nbC , nbT, nbG,nbN, nbGap)
				print( tab)
				exit()
			if nbGap == 0 :
				if nbO != sommeACTG and nbN == 0 :
					nbSNP+=1
					for lentabi in range(0,len(tab)):
						dicoSeqSNP[tableauSoucheName[lentabi]] += (tab[lentabi])

		nbSNPtotal=len(dicoSeqSNP[tableauSoucheName[lentabi]])
		if nbSNPtotal == 0:
			fileEmpty += 1
			listFileEmpty.append(filein)



					#print(lentabi)
				#dicoSeqSNP[]
			#print(tab)
				#print(nbO, sommeACTG)
		#print("il y a %i SNP" %nbSNP)
		#print(dict2fasta(dicoSeqSNP))
		#except:
			##print("Unexpected error:", sys.exc_info()[0])
			#fileEmpty += 1
			#listFileEmpty.append(filein)
			#print("ERREUR dans le fichier %s , pas d'alignement lisible (possible filtre 0 gap supprime toutes les positions)" % filein)
			#exit()





		listMotif = []

		for geneId, sequence in dicoSeqSNP.items():
			nbSNPtotal = (len(sequence))
			listMotif.append(sequence)

		nameMGG = filein.split("/")[-1].replace("_Orthologue_macse_NT.fasta","")
		if nameMGG not in dicoFileCountSNP.keys():
			dicoFileCountSNP[nameMGG] = {"NBSNP":nbSNPtotal,
										 "lenAlign":lenAlignement}


		dicoFilenbSNP[nbSNPtotal] = (dicoFilenbSNP.get(nbSNPtotal, 1))+1


		#print(nbSNPtotal)

		dicoCompteMotif = {k: listMotif.count(k) for k in set(listMotif)}
		#print(dict2txt(dicoCompteMotif))

		#test si 1 SNP et partager pas combien de souche:
		print("il y a %i SNP\n"% len(list(dicoCompteMotif)[0]))
		nbSNP = len(list(dicoCompteMotif)[0])

		if nbSNP == 2:
			values = list(dicoCompteMotif.values())
			if 1 in values:
				print("the file %s containe SNP Uninformatif\n" % filein)

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






	with open(outputfilename,"w") as output_handle:
		output_handle.write(dictList2txt(dicoOutputTxt))

	outputListEmpty = open("outputListEmpty.txt", "w")
	for fileEmptyName in listFileEmpty:
		outputListEmpty.write(fileEmptyName+"\n")
	outputListEmpty.close()

	with open("NBSNP-file.log","w") as output_handle:
		output_handle.write(dictDict2txt(dicoFileCountSNP))

	with open("NBSNP-file-count.log","w") as output_handle:
		output_handle.write(dict2txt(dicoFilenbSNP))

	print("\n\nExecution summary:")

	print("  - Outputting \n\
	Il y a au final %i SNP dans tout les MGG\n\
	Il y a %i fichiers vides\n\
	les sequences sont ajouter dans le fichier %s\n\
	la liste des fichiers vides est dans le fichier outputListEmpty.txt" %(nbSNPtotal,fileEmpty,outputfilename))
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
