#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package extractSeqFastaCorresponding.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath, existant_file, extractListFromFasta2, loadInDict2, directory

## Python modules
import argparse, os, subprocess
from time import localtime, strftime
import glob

## BIO Python modules
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
from MODULES_SEB import *

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='05/07/2016'
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
	parser = argparse.ArgumentParser(prog='extractSeqFastaCorresponding.py', description='''This Programme extract Fasta Seq with liste keep''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display extractSeqFastaCorresponding.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-f', '--fasta', metavar="<path/to/fasta>", type = directory, required=True,, dest = 'fastaFile', help = 'fasta files')
	filesReq.add_argument('-l', '--list', metavar="<path/to/files>", type = directory, required=True,v, dest = 'listFile', help = 'list files corresponding name')
	filesReq.add_argument('-o', '--out', metavar="<path/to/out>", type = directory, required=True,, dest = 'paramoutfile', help = 'Name of output file')
	filesReq.add_argument('-m', '--mgg', metavar="<filemane>", type=existant_file, required=True, dest = 'mggFileKeep', help = 'filename')

	# Check parameters
	args = parser.parse_args()
	checkParameters(args)

	#Welcome message
	print("#################################################################")
	print("#        Welcome in extractSeqFastaCorresponding (Version " + version + ")          #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	fastaFile = args.fastaFile
	listFile = args.listFile
	outputfilePath = args.paramoutfile
	mggFileKeep = relativeToAbsolutePath(args.mggFileKeep)


	print("\t - Path with fasta is: %s" % fastaFile.pathDirectory)
	print("\t - Path with corresponding Orthologues is  : %s" % listFile.pathDirectory)
	print("\t - MGG list keep are in file: %s\n" % mggFileKeep)

	print("\t - Output Orthologues fasta is: %s\n\n" % outputfilePath)






	if fastaFile[-1] != "/":
		fastaFile += "/"
	if outputfilePath[-1] != "/":
		outputfilePath += "/"
	if listFile[-1] != "/":
		listFile += "/"

	#os.system("rm "+outputfilePath+"*.fas")

	#recupération de la liste des CDS complet
	listCDSfiles = lsExtInDirToList(fastaFile,"fasta")

	#ouverture de la liste des MGG à garder
	mggKeepall = loadInList(mggFileKeep)
	# trie de la list pour supprimer les T1 et T2 avec T0
	MGGWithoutT = ["MGG_15255"]
	mggKeep = []
	count = 0
	toRM = []
	for mgg in mggKeepall:
		mggNoT = mgg.replace("T0","").replace("T1","").replace("T2","")
		if mggNoT not in MGGWithoutT:
			MGGWithoutT.append(mggNoT)
			mggKeep.append(mgg)
		else:
			count+=1
			toRM.append(mgg)

	with open("toRM.txt","w") as toRMFile:
		for rm in toRM:
			toRMFile.write(rm+"\n")
	with open("7037MGGothologue.txt","w") as mggKeepFile:
		for kepp in mggKeep:
			mggKeepFile.write(kepp+"\n")
	print("count = %i" % count)
	print("nbMGGkeep = %i" % len(mggKeep))

	 ##ouverture des fichiers de correspondances dans un dico de dico
	dicoCorrespondanceMGGTocontig = {}
	dicoCorrespondancecontigToMGG = {}
	listcorespFiles = lsExtInDirToList(listFile,"")
	#print(listcorespFiles)
	for file in listcorespFiles:
		name = file.split("/")[-1].split("_")[0]
		#print(name)
		dicoCorrespondanceMGGTocontig[name] = loadInDict(file)
		dicoCorrespondancecontigToMGG[name] = loadInDict2(file)

	# parcours des fichiers fasta
	for fileCDS in listCDSfiles:
		#try:
		fileName = fileCDS.split("/")[-1].split(".")[0]
		#print(fileName)
		# ouverture du fichier de sortie
		output_handle = open(outputfilePath+fileName+"_CDS_MGGorthologue.fas", "w")

		#construction de la liste à garder des mgg
		listFileKeep = []
		ctmp = 0
		for MGG in mggKeep:
			ctmp+=1
			listFileKeep.append(dicoCorrespondanceMGGTocontig[fileName][MGG][0])
		#print("listFileKeep = %i" % len(listFileKeep))
		#print("\n")
		#print(listFileKeep)
		# ouverture des sequence completes est extraction des sequences d'intérets
		dico_keep = extractListFromFasta2(fileCDS, listFileKeep)
		nbKeep = len(dico_keep.keys())
		print("il y a %i sequences dans le fichier %s" %(nbKeep,output_handle.name))


		for geneId, record in dico_keep.items():
			MGGName = dicoCorrespondancecontigToMGG[fileName][geneId.id]
			oldNumID = record.id
			new_record_name = MGGName+"_"+oldNumID
			record.id = new_record_name
			record.name = ""
			seq = record.seq
			SeqIO.write(record,output_handle, "fasta")
		#except:
			#print(fileName, MGG)
		#exit()
	output_handle.close()

	######
	os.system("rm "+outputfilePath+"orthologue/*.fasta")
	#Concatenation des orthologues de Farman et Gemo
	listFastaOut = lsFastaInDirToList(outputfilePath)
	for fastaFile in listFastaOut:
		dictSequences = fasta2dict(fastaFile)
		for geneId, record in dictSequences.items():
			#print(geneId)
			MGGName = "_".join(geneId.split("_")[0:2])
			if MGGName in toRM or "T1" in MGGName or "T2" in MGGName :
				#print(MGGName)
				MGGName = MGGName.replace("T0","T0").replace("T1","T0").replace("T2","T0")

			souche = geneId.split("_")[2]
			#print(souche)
		# ouverture du fichier de sortie
			output_handle = open(outputfilePath+"orthologue/"+MGGName+"_Orthologue.fasta", "a")
			new_record_name = souche
			record.id = ""
			record.id = new_record_name
			record.name = ""
			seq = record.seq
			SeqIO.write(record,output_handle, "fasta")
			output_handle.close()
	#listfastaDarren = lsFastaInDirToList("/media/sebastien/Bayer/hotespe/NO_FILTER/darren")
	#for fastaFile in listfastaDarren:
		#dictSequences = fasta2dict(fastaFile)
		#for geneId, record in dictSequences.items():
			#MGGName = geneId
			#if MGGName in mggKeep:
				#souche = fastaFile.split("/")[-1].split("_")[0]
			## ouverture du fichier de sortie
				#output_handle = open(outputfilePath+"orthologue/"+MGGName+"_Orthologue.fasta", "a")
				#new_record_name = souche
				#record.id = ""
				#record.id = new_record_name
				#record.name = ""
				#seq = record.seq
				#SeqIO.write(record,output_handle, "fasta")

	dico1,dico2 = nbSeqInFile2dict(outputfilePath+"orthologue/")
	print("check if NBsouche and sequences are correctly extract:\n")
	print(dict2txt(dico2))
	print("\n\nIf up are same below OK\n\n%i\t%i" % (count, len(mggKeep)))

	#print("  - Outputting \n\
	#Il y a au final %i Sequences garder\n\
	#les sequences sont ajouter dans le fichier %s" %(nbKeep,outputfilePath))
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
