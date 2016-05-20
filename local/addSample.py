#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
## @package addSample.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import directory, relativeToAbsolutePath, nbSeqInFile2dict, dict2txt,fasta2dict

## Python modules
import argparse
from time import localtime, strftime

## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='20/04/2016'
debug="False"
#debug="True"


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
# 	start=time.clock()
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='addSample.py', description='''This Programme find corresponding MGG and BR32 genes''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display addSample.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	files = parser.add_argument_group('Input info for running')
	files.add_argument('-s', '--samples', metavar="<path/to/filein>",type = directory, required=True, dest = 'pathSampleIn', help = 'path to fasta sample files in')
	files.add_argument('-c', '--csv', metavar="<filename>", dest = 'csvFileParam', help = 'Name of csv file with Orthologue infos')
	files.add_argument('-pi', '--pathin', metavar="<path/to/filein>",type = directory, required=True, dest = 'pathDirectoryIn', help = 'path to fasta files in')
	files.add_argument('-po', '--pathout', metavar="<path/to/fileout>",type = directory, required=True, dest = 'pathDirectoryOut', help = 'path to fasta files Out')
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#            Welcome in addSample (Version " + version + ")              #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	pathDirectoryIn = args.pathDirectoryIn
	pathDirectoryOut = args.pathDirectoryOut
	pathSampleIn = args.pathSampleIn

	dicoNbSeqInFiles,dicoNbFilesNbSouche = nbSeqInFile2dict(str(pathSampleIn.pathDirectory))
	print((dict2txt(dicoNbSeqInFiles)))
	print(dict2txt(dicoNbFilesNbSouche))

	# ajoute à la variable current_dir le chemin ou est executer le script
	current_dir = os.path.dirname(os.path.abspath(__file__))

	# Ouverture du fichier OrthoMCL pour correspondance famille nom gene
	dictNameFamilly={}
	with open(args.csvFileParam,"r") as fichier:
		for ligne in fichier:
			ltab = ligne.rstrip().split("\t")
			geneName=ltab[1]
			famillyName=ltab[2]
			if famillyName not in dictNameFamilly.keys():
				dictNameFamilly[famillyName]={	"MGG":"",
												"BR32":""
														}
			if "MGG" in geneName.split("_")[0]:
				dictNameFamilly[famillyName]["MGG"] = geneName

			elif "M_BR32" in "_".join(geneName.split("_")[0:2]):
				dictNameFamilly[famillyName]["BR32"] = geneName

		#print(dictNameFamilly)
		#for key in sorted(dictNameFamilly.keys()):
			#value = key+"\t"+str(dictNameFamilly[key])
			#print(value)


	#print(dictDict2txt(dictNameFamilly))
	dicoMGG2BR32 = {}
	for key, dico in dictNameFamilly.items():
		MGG = dico["MGG"]
		BR32 = dico["BR32"]

		if BR32 != "":
			#print("%s\t%s" % (MGG, BR32))
			dicoMGG2BR32[MGG]=BR32

	#print(dict2txt(dicoMGG2BR32))
	print(len(dicoMGG2BR32))

	# Récupère liste des fasta clean
	listFiles = ["_".join(filename.split("/")[-1].split("_")[0:2]) for filename in pathDirectoryIn.listFiles]
	#print(listFiles)

	notMGG = []
	BR32List = []
	for MGG in listFiles:
		if MGG in dicoMGG2BR32.keys():
			BR32List.append(dicoMGG2BR32[MGG])
		else:

			notMGG.append(MGG)

	print("BR32List count:%i" %(len(BR32List)))
	print("notMGG count:%i" %(len(notMGG)))


	dicoFastaSample = {}
	for sampleFile in pathSampleIn.listFiles:
		souche = sampleFile.split("/")[-1].split(".")[0]
		#print(souche)
		dicoFastaSample[souche] = fasta2dict(sampleFile)

	#print(dicoFastaSample)

	#print(pathDirectoryOut.pathDirectory)

	count = 0
	for filename in pathDirectoryIn.listFiles:
		basename = filename.split("/")[-1]
		MGG = "_".join(filename.split("/")[-1].split("_")[0:2])

		if MGG in dicoMGG2BR32.keys() :
			BR32ID = dicoMGG2BR32[MGG]

			first = 0
			for souche in dicoFastaSample.keys():
				#print(souche)
				if BR32ID in dicoFastaSample[souche].keys():

					#print(souche+"\t"+basename+"\t"+MGG+"\t"+dicoFastaSample[souche][BR32ID].id)

					if first == 0:
						count+=1
						try:
							outputFile.close()
						except:
							pass
						outputFile = open(pathDirectoryOut.pathDirectory+basename,"w")  # Open output file
						#	Ouverture des sequences fasta et chargement dans dictionnaire et ecriture d'un fichier out qui contiendra les sequences renommer
						record_dict = SeqIO.to_dict(SeqIO.parse(open(filename, "rU"), "fasta"))
						for ID, record in record_dict.items():
							SeqIO.write(record,outputFile, "fasta")
						first=1
					record = dicoFastaSample[souche][BR32ID]
					newRecord = SeqRecord(Seq(str(record.seq)),id=souche,name=souche, description=record.id)
					SeqIO.write(newRecord,outputFile, "fasta")
	try:
		outputFile.close()
	except:
		pass
	print(count)

	dicoNbSeqInFiles,dicoNbFilesNbSouche = nbSeqInFile2dict(str(pathDirectoryOut.pathDirectory))

	with open("dicoNbSeqInFiles.log","w") as output_handle:
		output_handle.write(dict2txt(dicoNbSeqInFiles))

	with open("dicoNbFilesNbSouche.log","w") as output_handle:
		output_handle.write(dict2txt(dicoNbFilesNbSouche))

	countCP = 0
	for fileName, nbSouche in dicoNbSeqInFiles.items():
		if int(nbSouche) == 48:
			os.system("cp "+pathDirectoryOut.pathDirectory+fileName+".fas "+pathDirectoryOut.pathDirectory+"/../test/")
			countCP+=1


	#alignDirecory = directory(current_dir+"/test/")
	#for nameFile in alignDirecory.listFiles:
		#clustalw_exe = "/usr/bin/clustalw"
		#clustalw_cline = ClustalwCommandline(clustalw_exe, infile=nameFile)
		#assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
		#stdout, stderr = clustalw_cline()




	print("\n\nExecution summary:")

	print("  - Outputting \n\
     il y a %i file copier" % countCP)

	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
