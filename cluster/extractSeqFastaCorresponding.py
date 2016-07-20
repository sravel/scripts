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
from MODULES_SEB import relativeToAbsolutePath, existant_file, extractListFromFasta, loadInDictCol, directory,dictDict2txt, nbSeqInFile2dict

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

def extractListFromFasta2(sequenceFile,FileList ):
	dicoOutput = {}
	# Ouverture des sequences fasta MGG et chargement dans dictionnaire
	dictSequences = fasta2dict(sequenceFile)

	# ouverture des identifiants a garder
	listKeep = FileList
	keep = 0
	noKeep = 0
	for ID, record in dictSequences.items():
		if ID in listKeep:
			keep +=1
			dicoOutput[record] = record
		else:
			noKeep += 1
	return dicoOutput

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
	parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-f', '--fasta', metavar="<path/to/fasta>", type = directory, required=True, dest = 'fastaFile', help = 'fasta files')
	filesReq.add_argument('-l', '--list', metavar="<path/to/files>", type = directory, required=True, dest = 'listFile', help = 'list files corresponding name')
	filesReq.add_argument('-o', '--out', metavar="<path/to/out>", type = directory, required=True, dest = 'paramoutfile', help = 'Name of output file')
	filesReq.add_argument('-m', '--mgg', metavar="<filemane>", type=existant_file, required=True, dest = 'mggFileKeep', help = 'filename')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#     Welcome in extractSeqFastaCorresponding (Version " + version + ")     #")
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




	#recupération de la liste des CDS complet
	listCDSfiles = fastaFile.lsExtInDirToList(["fasta", "fas", "fa"])
	print("\n".join(listCDSfiles))

	#ouverture de la liste des MGG à garder
	mggKeepall = loadInList(mggFileKeep)


	# trie de la list pour supprimer les T1 et T2 avec T0
	MGGWithoutT, mggKeep, toRM = [], [], []

	for mgg in mggKeepall:
		mggNoT = mgg.replace("T0","").replace("T1","").replace("T2","")
		if mggNoT not in MGGWithoutT:
			MGGWithoutT.append(mggNoT)
			mggKeep.append(mgg)
		else:
			toRM.append(mgg)

	with open("List"+str(len(toRM))+"TranscriptAlternatifstoRM.txt","w") as toRMFile:
		txt = "\n".join(toRM)
		toRMFile.write(txt)

	with open("List"+str(len(mggKeep))+"MGGothologuesKEEP.txt","w") as mggKeepFile:
		txt = "\n".join(mggKeep)
		mggKeepFile.write(txt)


	print("count = %i" % len(toRM))
	print("nbMGGkeep = %i" % len(mggKeep))



	 ##ouverture des fichiers de correspondances dans un dico de dico
	dicoCorrespondanceMGGTocontig = {}
	dicoCorrespondancecontigToMGG = {}

	listcorespFiles = listFile.lsExtInDirToList("")

	if args.debug == "True" : print("\n".join(listcorespFiles))

	for file in listcorespFiles:
		name = file.split("/")[-1].split("_")[0]
		#print(name)
		dicoCorrespondanceMGGTocontig[name] = loadInDictCol(file,0,1)
		dicoCorrespondancecontigToMGG[name] = loadInDictCol(file,1,0)

	# parcours des fichiers fasta
	#for fileCDS in listCDSfiles:

		#fileName = fileCDS.split("/")[-1].split(".")[0]

		## ouverture du fichier de sortie
		#with open(outputfilePath+fileName+"_CDS_MGGorthologue.fas", "w") as output_handle:

			##construction de la liste à garder des mgg
			#listFileKeep = []
			#ctmp = 0
			#for MGG in mggKeep:
				#ctmp+=1
				#listFileKeep.append(dicoCorrespondanceMGGTocontig[fileName][MGG])

			#dico_keep = extractListFromFasta2(fileCDS, listFileKeep)
			#nbKeep = len(dico_keep.keys())
			##print("il y a %i sequences dans le fichier %s" %(nbKeep,output_handle.name))


			#for geneId, record in dico_keep.items():
				#MGGName = dicoCorrespondancecontigToMGG[fileName][geneId.id]
				#oldNumID = record.id
				#new_record_name = MGGName+"_"+oldNumID
				#record.id = new_record_name
				#record.name = ""
				#seq = record.seq
				#SeqIO.write(record.upper(),output_handle, "fasta")


	######
	os.makedirs(outputfilePath+"orthologue/", exist_ok=True)
	os.system("rm "+outputfilePath+"orthologue/*.fasta")

	#Concatenation des orthologues de Farman et Gemo
	listFastaOut = lsFastaInDirToList(outputfilePath)
	nblignetotal = len(listFastaOut)
	print(listFastaOut)
	ctr = 0
	for fastaFile in listFastaOut:
		dictSequences = fasta2dict(fastaFile)
		print(fastaFile)
		percent = (float(ctr)/float(nblignetotal))*100
		sys.stdout.write("\rProcessed up to %0.2f %%..." % percent)
		sys.stdout.flush()

		for geneId, record in dictSequences.items():
			#print(geneId)
			MGGName = "_".join(geneId.split("_")[0:2])
			if MGGName in toRM or "T1" in MGGName or "T2" in MGGName :
				#print(MGGName)
				MGGName = MGGName.replace("T0","T0").replace("T1","T0").replace("T2","T0")

			souche = geneId.split("_")[2]
			#print(souche)
			# ouverture du fichier de sortie
			#with open(outputfilePath+"orthologue/"+MGGName+"_Orthologue.fasta", "a") as output_handle:
			output_handle = open(outputfilePath+"orthologue/"+MGGName+"_Orthologue.fasta", "a")
			new_record_name = souche
			record.id = ""
			record.id = new_record_name
			record.name = ""
			seq = record.seq
			SeqIO.write(record.upper(),output_handle, "fasta")
			output_handle.close()
		ctr+=1

	#dico1,dico2 = nbSeqInFile2dict(outputfilePath+"orthologue/")
	#print("check if NBsouche and sequences are correctly extract:\n")
	#print(dict2txt(dico2))
	#print("\n\nIf up are same below OK\n\n%i\t%i" % (count, len(mggKeep)))

	#print("  - Outputting \n\
	#Il y a au final %i Sequences garder\n\
	#les sequences sont ajouter dans le fichier %s" %(nbKeep,outputfilePath))
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
