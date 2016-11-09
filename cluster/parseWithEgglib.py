#!/usr/local/bioinfo/python/2.7.9_build2/bin/python
# -*- coding: utf-8 -*-
# @package parseWithEgglib.py
# @author Sebastien Ravel

"""
	The parseWithEgglib script
	==========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 13/07/2016
	:version: 0.1

	Script description
	------------------

	This Program takes fasta alignement and check SNP info\n
	Le script regarde avec Egglib la liste des sites polymorphes et les singletons:\n
	si la liste des sites polymorphes == 0 ou les 2 listes (singletons et sites polymorphes) sont identiques ajout dans le liste des CDS à remove

	Example
	-------

	>>> parseWithEgglib.py -f NT_ALIGN/ -o filterEgglib

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display parseWithEgglib version number and exit

	Input mandatory infos for running:
		- \-f <path/to/directory/fasta>, --fasta <path/to/directory/fasta>
						path to fasta files
		- \-o <path/to/directory>, --out <path/to/directory>
						Name of output file directory

	Input infos for running with default values:
		- \-i <filename>, --info <filename>
						output file with info results (default = sitesInformation.tab)
		- \-n <filename>, --noinfo <filename>
						output file list with gene with no info (default = outputListNoInfo.txt)

"""



##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import  dictDict2txt, directory

## Python modules
import argparse
from time import localtime, strftime
from egglib import Align

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='13/07/2016'

def compare_intersect(x, y):
	return list(set(x).intersection(y))

##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='parseWithEgglib.py', description='''This Program takes fasta check if SNP info''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
											'display parseWithEgglib version number and exit')

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-f', '--fasta', metavar="<path/to/directory/fasta>", type = directory, required=True, dest = 'fastaFileDir', help = 'path to fasta files')
	filesReq.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'pathOut', help = 'Name of output file directory')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-i', '--info', metavar="<filename>",type = str, default="sitesInformation.tab", required=False, dest = 'infoFile', help = 'output file with info results (default = sitesInformation.tab)')
	files.add_argument('-n', '--noinfo', metavar="<filename>",type = str, default="outputListNoInfo.txt", required=False, dest = 'noInfoFile', help = 'output file list with gene with no info (default = outputListNoInfo.txt)')


	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#           Welcome in parseWithEgglib (Version " + version + ")            #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')


	# Récupère les arguments
	pathFastaFile = args.fastaFileDir
	pathFileOut = args.pathOut

	# defaults option
	infoFile=pathFileOut.pathDirectory+args.infoFile
	noInfoFile=pathFileOut.pathDirectory+args.noInfoFile

	# resume value to user
	print(" - Intput Info:")
	print("\t - Working in directory: %s" % pathFileOut.pathDirectory)
	print("\t - Fasta were in directory: %s" % pathFastaFile.pathDirectory)

	print(" - Output Info:")
	print("\t - Output file with sites information: %s" % infoFile)
	print("\t - Output file list with genes no info: %s" % noInfoFile)


	dicoOutput = {}
	listFasta = pathFastaFile.lsExtInDirToList(["fasta","fas","fa","fna"])
	listFasta = [string.encode("utf-8") for string in listFasta]
	#print listFasta
	for filein in listFasta:
		name = filein.split("/")[-1]

		#namefile = filein.split("/")[-1]
		#align = Align(filein)
		#listIndexRM = []
		##print len(align.names())
		#count=0
		#for name in align.names():
			#nameS =  name.split(" ")[0]
			#if  "T86" not in nameS and "Mycfi" not in nameS:
				#count+=1
				#index = align.find(name)
				#del align[index]
		#x = Align.polymorphism(align)
		##print dict2txt(x)
		#dicoOutput[namefile] = {	"singletons":x['singletons'],
								#"siteIndices": x['siteIndices'],
								#"nbSNP": len(x['siteIndices'])}

		#print(name)
		align = Align(filein)

		x = Align.polymorphism(align)

		#print dict2txt(x)

		dicoOutput[name] = {	"singletons":x['singletons'],
								"siteIndices": x['siteIndices'],
								}


	with open(infoFile,"w") as output_handle:
		output_handle.write(dictDict2txt(dicoOutput))


	listFileRM = []
	for fileName, dico in dicoOutput.items():
		listSingleton = list(dico["singletons"])
		listSiteIndices = list(dico["siteIndices"])
		intersection = compare_intersect(listSingleton,listSiteIndices)
		if len(listSiteIndices) == 0 :
			listFileRM.append(fileName)
		elif len(compare_intersect(listSingleton,listSiteIndices))==len(listSiteIndices):
			listFileRM.append(fileName)
			#print listSingleton, listSiteIndices, len(listSingleton), len(listSiteIndices)





	with open(noInfoFile, "w") as outputListNoInfo:
		txt = "\n".join(listFileRM)
		outputListNoInfo.write(txt)



	#print dictDict2txt(dicoOutput)




	#print("\n\nExecution summary:")

	#print("  - Outputting \n\
	#Il y a au final %i Sequences gardées sur les %i initial\n\
	#les sequences sont ajoutées dans le fichier %s" %(nbKeep,nbTotal,outputfilename))
	#print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("\n\n#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
