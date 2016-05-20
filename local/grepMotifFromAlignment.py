#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
## @package grepMotifFromAlignment.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
#sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import directory, dictList2txt, dictDict2txt, dict2txt

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
def checkParameters (arg_list):
	# Check input related options
	if (not arg_list.pathDirectory):
		#print ('Error: No input file defined via option -i/--input !' + "\n")
		parser.print_help()
		parser.print_usage()
		parser.format_usage()
		parser.format_help()
		exit()


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='grepMotifFromAlignment.py', description='''This Programme parse Aligment info''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display grepMotifFromAlignment.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	files = parser.add_argument_group('Input info for running')
	files.add_argument('-d', '--directory', metavar="<path/to/directory>", dest = 'pathDirectory', help = 'path to directory fasta files')
	files.add_argument('-o', '--out', metavar="<filename>", dest = 'paramoutfile', help = 'Name of output file')

	# Check parameters
	args = parser.parse_args()
	checkParameters(args)

	#Welcome message
	print("#################################################################")
	print("#        Welcome in grepMotifFromAlignment (Version " + version + ")          #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	pathDirectory = directory(args.pathDirectory)

	outputfilename = args.paramoutfile

	# ouverture de la taille des MGG dans dico
	#dicoLenMGG = loadInDict("len_6878MGG.txt")

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
			if int(allcount) != 46:
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
