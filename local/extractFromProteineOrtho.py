#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
## @package extractFromProteineOrtho.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
## Python modules
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import directory, existant_file, dict2txt, sort_human, dictDict2txt, AutoVivification,printCol

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
VERSION_DATE='06/06/2016'
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
	parser = argparse.ArgumentParser(prog='extractFromProteineOrtho.py', description='''This Programme take proteineOrtho output and take file with Orthologue 1/1''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display extractFromProteineOrtho version number and exit')
	parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	files = parser.add_argument_group('Input info for running')
	files.add_argument('-po', '--pathout', metavar="<path/to/fileout>", type = directory, required=True, dest = 'pathFileOut', help = 'path to fasta files Out')
	files.add_argument('-p', '--proteine', metavar="<filename>",type=existant_file, required=True, dest = 'proteineOrthoFile', help = 'proteineOrthoFile')
	files.add_argument('-r', '--ref', metavar="<string>", required=True, dest = 'refName', help = 'Name of strain reference (ex: Mycfi) ')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#      Welcome in extractFromProteineOrtho (Version " + version + ")        #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument

	pathFileOut = args.pathFileOut
	ref = args.refName
	workingDir = "/".join(pathFileOut.pathDirectory.split("/")[:-2])+"/"
	correspondingCDSDir = workingDir+"correspondingCDS-contig/"

	print("\t - Input Path is: %s" % pathFileOut.pathDirectory)
	print("\t - ref strain is : %s" % ref)
	print("\t - Working directory is: %s" % workingDir)
	print("\t - Corresonding CDS ref/strain directory is: %s\n\n" % correspondingCDSDir)

	# liste de toute les souches de proteineOrtho
	listSouches =[]

	# dico de proteine orthologue
	dico_ortho = {}
	exist=0
	# creer le répertoire contenant les correspondance entre ref et souches
	if os.path.exists(correspondingCDSDir):
		printCol.yellow("  Warning , folder "+correspondingCDSDir+" already exist !!!!" )
		exist=1
	else:
		os.makedirs(correspondingCDSDir)
	if exist == 1:
		printCol.yellow("  Do you want to remove all analysis? (y/n)\n")
		inp = None
		while inp == None and inp != "y" and inp != "n" and inp != "yes" and inp != "no":
			inp = input()
			if inp == "y":
				os.popen('rm -r '+correspondingCDSDir+'; mkdir '+correspondingCDSDir )
				exist=0

			elif inp == "n":
				printCol.red(">>>Program exit\n")
				exit()


	# ouverture du fichier de résultat protine ortho pour construire liste de seq ortho
	with open(args.proteineOrthoFile,"r") as proteineOrthoFile:
		for line in proteineOrthoFile:
			if "#" not in line:
				tabline = line.split("\t")
				souche_contig1 = tabline[0]						# nom de la souche et contig
				souche_contig2 = tabline[1]
				namesouche1 = tabline[0].split("_")[0]		# nom de la souche seul
				namesouche2 = tabline[1].split("_")[0]

				if args.debug == "True" : print(souche_contig1, namesouche1,souche_contig2, namesouche2)				# DEBUG

				if ref in namesouche1:
					dico_ortho.setdefault(souche_contig1, []).append(souche_contig2)
					if namesouche2 not in listSouches:
						listSouches.append(namesouche2)

				if ref in namesouche2:
					dico_ortho.setdefault(souche_contig2, []).append(souche_contig1)
					if namesouche1 not in listSouches:
						listSouches.append(namesouche1)


	listSouchessort = sorted(listSouches, key=sort_human)
	if args.debug == "True" : print(dict2txt(dico_ortho))
	if args.debug == "True" : print(listSouchessort)
	if args.debug == "True" : print(len(listSouchessort))

	# ecriture des correspondence othologue 1/1:
	nbOrtho1_1=0
	for souche_contig1 , listcorresp in dico_ortho.items():
		namesouche1 = souche_contig1.split("_")[0]
		tabsoucheFind = [souche.split("_")[0] for souche in listcorresp]

		if len(listcorresp) == len(listSouchessort) and sorted(listSouchessort, key=sort_human) == sorted(tabsoucheFind, key=sort_human):
			for soucheFind in listcorresp:
				souche_contig2 = soucheFind
				namesouche2 = souche_contig2.split("_")[0]
				correspondanceMGGContig = open(correspondingCDSDir+namesouche2+"_corespondingMGG-contig","a")
				correspondanceMGGContig.write("%s\t%s\n"%(souche_contig1,souche_contig2))
			nbOrtho1_1+=1

	print("\t - %i orthologues 1/1 found on the %s strains" % (nbOrtho1_1, len(listSouchessort)+1))

	# ouverture d'un tableau résumer
	with open(workingDir+"Orthologue_Tab_Stats.tab","w") as tabFileOut:

		dicoCountNB = AutoVivification()
		for souche_contig1 in sorted(dico_ortho.keys(), key=sort_human):
			namesouche1 = souche_contig1.split("_")[0]

			# initialisation du dico a 0
			for souche in listSouchessort:
				dicoCountNB[souche_contig1][souche] = 0

			for souche_contig2 in dico_ortho[souche_contig1]:
				namesouche2 = souche_contig2.split("_")[0]

				dicoCountNB[souche_contig1][namesouche2]+=1


		tabFileOut.write(dictDict2txt(dicoCountNB,ref))






	print("\n\nExecution summary:")

	print("  - Outputting \n\
     - %s\n\
     - %s\n\n" % (tabFileOut.name,correspondingCDSDir) )

	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")












	#for fileCDS in listCDSfiles:
		#fileName = fileCDS.split("/")[-1].split("_")[0]
		#listFile.write(current_dir+pathFileOut+fileName+".fasta\n")

		#print(fileName)
		#if nbstart !=0 or nbstartandstop !=0 or nbstop !=0:
			#total = nbstart+nbstop+nbstartandstop
			#statFile.write("%s\t%s\t%s\t%s\t%s\n" %(fileName,nbstart,nbstop, nbstartandstop,total))
			#nbstart=0
			#nbstop=0
			#nbstartandstop=0

		#output_file = open(pathFileOut+fileName+".fasta", "w")

		#record_dict = fasta2dict(fileCDS)
		#for name in sorted(record_dict.keys()):
			#record = record_dict[name]
			#oldNumID = record.id
			#new_record_name = fileName+"_"+oldNumID
			#record.id = new_record_name
			#record.name = ""
			#seq = record.seq
			#firstCodon = seq[:3]
			#endCodon = seq[-3:]

			## Test ATG start and Stop codons
			#if str(firstCodon.upper()) in "ATG":
				#startATG=1

			#else:
				#startATG=0
			#if str(endCodon.upper()) in ["TAG","TAA","TGA"]:
				#stop=1
			#else:
				#stop=0
			#if startATG == 1 and stop == 1:
				#nbstartandstop+=1
				#record.description = 'start_stop'
			#elif startATG == 1 and stop == 0:
				#nbstart+=1
				#record.description = 'start_-'
			#elif startATG == 0 and stop == 1:
				#nbstop+=1
				#record.description = '-_stop'


			## write the whole thing out
			#SeqIO.write(record, output_file, 'fasta')
