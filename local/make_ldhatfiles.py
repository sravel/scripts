#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
## @package make_ldhatfiles.py
# @author Lea Picard

##################################################
## Modules
##################################################
## Python modules
import argparse, re, os
from subprocess import check_output
import egglib # USE EGGLIB_3

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='15/03/2016'

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

		if nameChro not in dictListPositions.keys():
			dictListPositions[nameChro] = [posSNP]
		else:
			dictListPositions[nameChro].append(posSNP)

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
	outputfilename = paramfilename.split(".")[0]+".sites"
	outputfile = open(outputfilename, 'w')

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

	outputfile.write("%i %i %i\n" %(nbInd, nbSNP, dataType))

	for souche in orderlist:
		outputfile.write(">"+souche+"\n"+dictseq[souche]+"\n")

	outputfile.close()

	####### Remove temps file
	os.remove(genometab)
	os.remove("temps1.tab")
	return nbSNP, dictListPositions[nameChro], outputfile.name

def relativeToAbsolutePath(relative):
	from subprocess import check_output
	if relative[0] != "/":			# The provided path is a relative path, ie does not start with /
		command = "readlink -m "+relative
		absolutePath = check_output(command, shell=True).decode("utf-8").rstrip()
		return absolutePath+"/"
	else:						# Relative is in fact an absolute path, send a warning
		absolutePath = relative;
		return absolutePath+"/"


def dictDict2txt(dico):
	"""
	Function that takes a dictionary and returns a tabular string with::

		"key\\tvalue\\n".

	:param dico: a python dictionary
	:type dico: dict()
	:rtype: str()
	:return: string with "key\\tvalue\\n

	Example:
		>>> dico = {"Souche1":{"NUM":"171","MIN":"2042","MAX":"3133578","N50 BP":"938544","N50 NUM":"11"},
				    "Souche2":{"NUM":"182","MIN":"5004","MAX":"74254","N50 BP":"45245","N50 NUM":"45"}}
		>>> dictDict2txt(dico)
		Info	NUM	MIN	MAX	N50 BP	N50 NUM
		Souche1	171	2042	3133578	938544	11
		Souche2	182	5004	74254	45245	45

	"""

	txtoutput = ""
	headerc=0
	for key in sorted(dico.keys()):
		dicoInfosValues = [str(item) for item in dico[key].values()]
		if headerc == 0:
			header = dico[key].keys()
			value = "Info\t" + "\t".join(header)
			txtoutput += "%s\n" % str(value)
			headerc=1

		value = "\t".join(dicoInfosValues)
		txtoutput += "%s\t%s\n" % (str(key),str(value))
	return txtoutput

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
	files = parser.add_argument_group('Input info for running')
	files.add_argument('-wd', '--workdir', metavar="<path>", required=True, dest = 'workdir', help = 'Path of the directory where files will be created')
	files.add_argument('-t', '--tab', metavar="<filename>", required=True, dest = 'tabFile', help = 'Name of tab file in (input whole path if file is not in the current working directory')
	files.add_argument('-st', '--size_tab', metavar="<filename>", required=True, dest = 'sizeTab', help = 'Name of a tab file containing the identifiers of the subunits of division (chromosome/scaffold/contig) and their total size. If some scaffolds are not wanted, comment the line.')
	files.add_argument('-dt', '--datatype', metavar="<int>", default=1, type=int, choices=[1,2], dest = 'datatype', help = '1 for haplotypic data (default), 2 for genotypic')
	files.add_argument('-f', '--flag', metavar="<char>", default="L", choices=["L","C"], dest = 'flag', help = 'L for CO (default), C pour gene conversion')

	# check parameters
	args = parser.parse_args()

	# get arguments
	workdir = relativeToAbsolutePath(args.workdir)
	tabFile = args.tabFile
	sizeTab = args.sizeTab
	dataType = args.datatype
	flag = args.flag

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
				subdir = workdir+basename+"/"+chro
				if not os.path.exists(subdir):
					os.makedirs(subdir)

				outputName = workdir+basename+"/"+chro+"/"+basename+"_"+chro+".tab"

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
		dictNbSNP[chroName], listPos, fasta = build_sites(name, dataType)
		listFasta.append(fasta)
		dictListPos[chroName] = listPos

	# for each subunit and its list of SNP positions
	for checkChro, listPos in dictListPos.items():
		if checkChro in dictFilesOut.keys():
			outputLocsName = workdir+basename+"/"+checkChro+"/"+basename+"_"+checkChro+".locs"

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
	cs.add_stat('Pi')
	cs.add_stat('thetaW')


	# load alignement
	for nameFasta in listFasta:

		scaffold = nameFasta.split("/")[-1].split(".")[0].replace(basename+"_","")

		# use egglib
		align = egglib.io.from_fasta(nameFasta)
		stats = cs.process_align(align)		# extract polymorphism data

		# get number of SNPs in file
		nbSNP = stats['ls']


		# print results
		if scaffold not in dictThetaInfo:
			dictThetaInfo[scaffold] = {	"Theta_SNP":stats['thetaW'],
										"Pi":stats['Pi'],
										"Nb_SNPs":nbSNP,
										"Theta_allSNPs":stats['thetaW']*nbSNP,
										"Theta_scaffold":stats['thetaW']*nbSNP/int(dictSizes[scaffold])
										}

	dicoMeanTheta = {}
	sommeTheta,sommeSize  = 0, 0
	for scaffold, dico in dictThetaInfo.iteritems():
		sommeTheta += dico["Theta_allSNPs"]
		sommeSize += int(dictSizes[scaffold])

	thetaCoreGenome = sommeTheta/sommeSize

	with open(workdir+basename+"_ThetaValues.tab", "w") as ThetaTab:
		ThetaTab.write(dictDict2txt(dictThetaInfo))
		ThetaTab.write("\nthetaCoreGenome\t%f.5" % thetaCoreGenome)
