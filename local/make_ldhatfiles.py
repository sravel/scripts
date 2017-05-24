#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
# @package make_ldhatfiles.py
# @author Lea Picard, Sebastien RAVEL

"""
	The make_ldhatfiles script
	==========================
	:author: Sebastien Ravel, Lea Picard
	:contact: sebastien.ravel@cirad.fr
	:date: 08/07/2016
	:version: 0.1

	Script description
	------------------

	This Program takes a tab file and returns LDhat .sites and .locs files

	Example
	-------

	>>> make_ldhatfiles.py -wd outPath -t SNP_table.tab -st chomosomeSize.txt

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display make_ldhatfiles.py version number and exit

	Input mandatory infos for running:
		- \-wd <path>, --workdir <path>
						Path of the directory where files will be created
		- \-t <filename>, --tab <filename>
						Name of tab file in (input whole path if file is not
						in the current working directory
		- \-st <filename>, --size_tab <filename>
						Name of a tab file containing the identifiers of the
						subunits of division (chromosome/scaffold/contig) and
						their total size. If some scaffolds are not wanted,
						comment the line.

	Input infos for running with default values:
		- \-dt <int>, --datatype <int>
						1 for haplotypic data (default), 2 for genotypic
		- \-m <char>, --methode <char>
						rhomap or interval (default)
		- \-f <char>, --flag <char>
						L for CO (default), C pour gene conversion

"""

##################################################
## Modules
##################################################

## Python modules
from sys import version_info, version
try:
	assert version_info <= (3,0)
except AssertionError:
	print("You are using version %s but version 2.7.x is require for this script!\n" % version.split(" ")[0])
	exit(1)


#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import directory, relativeToAbsolutePath, dictDict2txt, existant_file

import argparse


try:
	import egglib3 as egglib # USE EGGLIB_3
	if int(egglib.version.split(".")[0]) != 3 :
		print("You are using not use egglib V3!\n" )
		exit(1)
except ImportError:
	print("You are not able to load egglib V3!\n" )

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Variables Globales
version = "0.1"
VERSION_DATE = '15/03/2016'
completeLDhatPATH = "completeLDhat"
#intervalLDhatPATH = "interval"
##intervalLDhatPATH = "rhomap"
statLDhatPATH = "statLDhat"


##################################################
## Functions
def build_sites(paramfilename, dataType):
	"""fonction adaptée du script build_concensusV2.py pour fichier .sites LDHAT"""

	fichier = open(paramfilename,"r")
	outfile = open("temps1.tab", "w")

	head = fichier.readline()
	outfile.write(head.replace("CHROM\t",""))

	dictListPositions = {}

	for ligne in fichier:
		ligne = ligne.rstrip()

		lligne = ligne.split("\t")

		nameChro = lligne[0]
		posSNP = lligne[1]

		dictListPositions.setdefault(nameChro, [posSNP]).append(posSNP)

		#if nameChro not in dictListPositions.keys():
			#dictListPositions[nameChro] = [posSNP]
		#else:
			#dictListPositions[nameChro].append(posSNP)

		ligneoutput="\t".join(lligne[1:2])
		ref=lligne[2]

		souchestr="\t".join(lligne[3:]).replace("R",ref)
		ligneoutput+="\t"+lligne[2]+"\t"+souchestr+"\n"
		outfile.write(ligneoutput)
	outfile.close()
############################################################################################
#  switch matrice
############################################################################################
	fichier = open("temps1.tab","r")
	outfile = open("temps.tab", "w")
	A = []

	for ligne in fichier:
		tabligne = ligne.rstrip().split("\t")
		A.append(tabligne)
	#print(A)

	for ligne in list(zip(*A)):
		outfile.write("\t".join(ligne)+"\n")
	outfile.close()
############################################################################################
#  Grep consensus
############################################################################################
	# Récupère le fichier de conf passer en argument
	genometab = "temps.tab"

	if paramfilename.count(".") > 1:
		outputfilenameSite = ".".join(paramfilename.split(".")[0:-1])+".sites"
		outputfilenameFasta = ".".join(paramfilename.split(".")[0:-1])+".fasta"
	else:
		outputfilenameSite = paramfilename.split(".")[0]+".sites"
		outputfilenameFasta = paramfilename.split(".")[0]+".fasta"

	outputfileSite = open(outputfilenameSite, 'w')
	outputfileFasta = open(outputfilenameFasta, 'w')

	# Utilisation du tab
	TabFile = open(genometab, "r")
	dictseq={}
	head=TabFile.readline()

	orderlist=[]
	for tabline in TabFile:
		ltab=tabline.rstrip().split("\t")
		souche=ltab[0]
		if souche not in dictseq.keys():
			dictseq[souche]=""
			orderlist.append(souche)
			#get nb of sequences to add to file header
			nbInd = len(orderlist)

			seqreftab=ltab[1:]
			dictseq[souche]="".join(seqreftab)
			#get nb of SNPs in fasta sequence
			nbSNP = len(dictseq[souche])

	outputfileSite.write("%i %i %i\n" %(nbInd, nbSNP, dataType))

	for souche in orderlist:
		IDname = souche
		seq = dictseq[souche]
		record = SeqRecord(Seq(seq),id=IDname,name=IDname, description="")
		SeqIO.write(record,outputfileSite, "fasta")
		SeqIO.write(record,outputfileFasta, "fasta")


	outputfileSite.close()
	outputfileFasta.close()

	####### Remove temps file
	os.remove(genometab)
	os.remove("temps1.tab")
	return nbSNP, dictListPositions[nameChro], str(outputfileFasta.name), nbInd


############
## Main code
############
if __name__ == "__main__":

##
# parameters recovery
##

	parser = argparse.ArgumentParser(prog='make_ldhatfiles.py', description='''This Program takes a tab file and returns LDhat .sites and .locs files''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display make_ldhatfiles version number and exit')
	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-wd', '--workdir', metavar="<path>",type=directory, required=True, dest = 'workdir', help = 'Path of the directory where files will be created')
	filesreq.add_argument('-t', '--tab', metavar="<filename>",type=existant_file, required=True, dest = 'tabFile', help = 'Name of tab file in (input whole path if file is not in the current working directory')
	filesreq.add_argument('-st', '--size_tab', metavar="<filename>",type=existant_file, required=True, dest = 'sizeTab', help = 'Name of a tab file containing the identifiers of the subunits of division (chromosome/scaffold/contig) and their total size. If some scaffolds are not wanted, comment the line.')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-dt', '--datatype', metavar="<int>", default=1, type=int, choices=[1,2], dest = 'datatype', help = '1 for haplotypic data (default), 2 for genotypic')
	files.add_argument('-m', '--methode', metavar="<char>", default="interval", choices=["interval","rhomap"], dest = 'methode', help = 'rhomap or interval (default)')
	files.add_argument('-f', '--flag', metavar="<char>", default="L", choices=["L","C"], dest = 'flag', help = 'L for CO (default), C pour gene conversion')

	# check parameters
	args = parser.parse_args()


	#Welcome message
	print("#################################################################")
	print("#          Welcome in make_ldhatfiles (Version " + version + ")            #")
	print("#################################################################")
	# get arguments
	workingObjDir = args.workdir
	tabFile = relativeToAbsolutePath(args.tabFile)
	sizeTab = relativeToAbsolutePath(args.sizeTab)
	dataType = args.datatype
	intervalLDhatPATH = args.methode
	flag = args.flag



	print("\t - Workink Directory: %s" % workingObjDir.pathDirectory)
	print("\t - Input Path matrice is: %s" % tabFile)
	print("\t - Input Path size is: %s" % sizeTab)
	print("\t - dataType is : %s" % dataType)
	print("\t - Working with : %s" % intervalLDhatPATH)
	print("\t - flag is: %s\n\n" % flag)

	#exit()
##
# code
##

	# get basename to build the rest of the filenames
	basename = tabFile.split("/")[-1].split(".")[0]

	# build dictionary of scaffolds for the file to be split into
	dictSizes = {}

	with open(sizeTab, "r") as sizeTabOpen:
		for sizeLine in sizeTabOpen:
			# keys = IDs - values = total size
			checkChro = sizeLine.split("\t")[0]
			sizeChro = sizeLine.rstrip().split("\t")[1]

			dictSizes[checkChro] = sizeChro

	listRange = dictSizes.keys()

	## split by specified subunits (scaffold/chromosome/contig etc)

	# keys = subunits to be split, values = files to be written in
	dictFilesOut = {}

	with open(tabFile, "r") as tabOpen:
		# get the header of the original tab file to rebuild split tab file
		header = tabOpen.readline()

		# start from second line
		for line in tabOpen:
			# chro = identifier of the subunit in the first column
			chro = line.rstrip().split("\t")[0]

			# if chro considered belongs to user-defined range
			if chro in listRange:

				# create subdirectory for the current scaffold
				subdir = workingObjDir.pathDirectory+basename+"/"+chro
				if not os.path.exists(subdir):
					os.makedirs(subdir)

				outputName = workingObjDir.pathDirectory+basename+"/"+chro+"/"+basename+"_"+chro+".tab"

				# if chro not encountered yet, create file add header and current line
				if chro not in dictFilesOut.keys():
					dictFilesOut[chro] = open(outputName, "w")
					dictFilesOut[chro].write(header)
					dictFilesOut[chro].write(line)

				# otherwise just add current line to relevant file
				else:
					dictFilesOut[chro].write(line)

	# keys = names of split files, values = nb of SNPs in said file
	dictNbSNP = {}
	dictListPos = {}
	listFasta = []

	# for each split file
	for fileOut in dictFilesOut.values():
		name = fileOut.name
		chroName = name.split("/")[-1].split(".")[0].replace(basename+"_","")
		fileOut.close()
		# create corresponding .sites file and associate Nb of SNPs
		dictNbSNP[chroName], listPos, fasta, nbInd = build_sites(name, dataType)
		listFasta.append(fasta)
		dictListPos[chroName] = listPos

	# for each subunit and its list of SNP positions
	for checkChro, listPos in dictListPos.items():
		if checkChro in dictFilesOut.keys():
			outputLocsName = workingObjDir.pathDirectory+basename+"/"+checkChro+"/"+basename+"_"+checkChro+".locs"

			# create .locs file
			outputLocs = open(outputLocsName, "w")
			# write header as NbSNP ScaffSize Flag
			txt = "%i %s %s\n" %(dictNbSNP[checkChro], dictSizes[checkChro], flag)
			outputLocs.write(txt)

			# write SNP positions underneath
			txtLocs = " ".join(dictListPos[checkChro])+"\n"
			outputLocs.write(txtLocs)

			outputLocs.close()

	## calculate Pi and Theta values
	dictThetaInfo = {}

	cs = egglib.stats.ComputeStats()
	cs.add_stats('Pi','thetaW','ls_o')

	# load alignement
	for nameFasta in listFasta:

		scaffold = nameFasta.split("/")[-1].split(".")[0].replace(basename+"_","")

		# use egglib
		align = egglib.io.from_fasta(nameFasta, groups=False)
		stats = cs.process_align(align)		# extract polymorphism data

		# get number of SNPs in file
		nbSNP = align.ls
		#nbSNP = stats['ls_o']

		# print results
		if scaffold not in dictThetaInfo:
			dictThetaInfo[scaffold] = {	"Theta_SNP":stats['thetaW']/align.ls,
										"Pi":stats['Pi']/align.ls,
										"Nb_SNPs":nbSNP,
										"Theta_allSNPs":stats['thetaW'],
										"Theta_scaffold":stats['thetaW']/int(dictSizes[scaffold])
										}

	dicoMeanTheta = {}
	sommeTheta,sommeSize  = 0, 0
	for scaffold, dico in dictThetaInfo.iteritems():
		sommeTheta += dico["Theta_allSNPs"]
		sommeSize += int(dictSizes[scaffold])

	thetaCoreGenome = sommeTheta/sommeSize

	with open(workingObjDir.pathDirectory+basename+"/"+basename+"_ThetaValues.tab", "w") as ThetaTab:
		ThetaTab.write(dictDict2txt(dictThetaInfo))
		ThetaTab.write("\nthetaCoreGenome\t%.4f" % thetaCoreGenome)



	#MAKE sh script to run LDhat
	objDir = directory(workingObjDir.pathDirectory+basename)		# list all directory and files in the path


	#nbInd = 13
	#thetaCoreGenome = 0.007

	cmdLoadR = "module load compiler/gcc/4.9.2 bioinfo/geos/3.4.2 bioinfo/gdal/1.9.2 mpi/openmpi/1.6.5 bioinfo/R/3.2.2"
	cmdLookTable = completeLDhatPATH+" -n "+str(nbInd)+" -rhomax 100 -n_pts 201 -theta "+str(thetaCoreGenome)+" -prefix "+objDir.pathDirectory+basename

	with open(workingObjDir.pathDirectory+basename+"/runLDhat_"+basename+".sh", "w") as runSHFile:

		runSHFile.write("%s\n" % cmdLoadR)
		runSHFile.write("%s\n" % cmdLookTable)
		for scaff in sorted(objDir.listDir):
			scaffObjDir = directory(scaff)
			#print scaffObjDir.__repr__
			siteFile =	 [s for s in scaffObjDir.listFiles if ".site" in s][0]
			locsFile =	 [s for s in scaffObjDir.listFiles if ".locs" in s][0]
			basenameScaff = siteFile.split("/")[-1].split(".")[0]
			#print basename


			cmdCD = "cd "+scaff

			if "rhomap" in intervalLDhatPATH:
				cmdInterval = intervalLDhatPATH+" -seq "+siteFile+" -loc "+locsFile+" -lk "+objDir.pathDirectory+basename+"new_lk.txt -its 5000000 -bpen 10 -burn 100000 -samp 5000 -prefix "+scaffObjDir.pathDirectory+basenameScaff

			if "interval" in intervalLDhatPATH:
				cmdInterval = intervalLDhatPATH+" -seq "+siteFile+" -loc "+locsFile+" -lk "+objDir.pathDirectory+basename+"new_lk.txt -its 5000000 -bpen 10 -samp 5000 -prefix "+scaffObjDir.pathDirectory+basenameScaff

			cmdStat = statLDhatPATH+" -input "+scaffObjDir.pathDirectory+basenameScaff+"rates.txt -prefix "+scaffObjDir.pathDirectory+basenameScaff


			cmdGraph = "makeLDhatgraphs.R -f "+scaffObjDir.pathDirectory+basenameScaff+"rates.txt -o "+scaffObjDir.pathDirectory+basenameScaff+""




			#print "%s\n%s\n%s\n%s\n" % (cmdCD,cmdLookTable,cmdInterval,cmdStat)

			runSHFile.write("%s\n%s\n%s\n%s\n" % (cmdCD,cmdInterval,cmdStat,cmdGraph))





	os.system("chmod 755 "+workingObjDir.pathDirectory+basename+"/runLDhat_"+basename+".sh")


	cmdQsub = "qsub -V -q long.q -N "+basename+" -b Y -pe parallel_smp 4 "+workingObjDir.pathDirectory+basename+"/runLDhat_"+basename+".sh"

	print(cmdQsub)


	#print("\n\nExecution summary:")

	#print("  - Outputting \n\
	#- %s\n\
	#- %s\n\
	#- %s\n\n" % (tabFileOut.name,listKeepFile.name,correspondingCDSDir) )


	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
