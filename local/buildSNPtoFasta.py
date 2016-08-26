#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
# @package buildSNPtoFasta.py
# @author Sebastien Ravel

"""
	The buildSNPtoFasta script
	==========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 08/07/2016
	:version: 0.1

	Script description
	------------------

	This Programme parse GFF, TAB to add SNP in sequences

	Example
	-------

	>>> buildSNPtoFasta.py -g Myfi.gff3 -l List3873MGGothologuesKEEP.txt -t 62souches_nofilter.tab -f orthologue -o out

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display buildSNPtoFasta.py version number and exit

	Input mandatory infos for running:
		- \-g <filename>, --gff <filename>
						gff files with annotation
		- \-l <filename>, --list <filename>
						File with geneID to keep
		- \-t <filename>, --tab <filename>
						File with SNP
		- \-f <path/to/directory>, --fasta <path/to/directory>
						Directory with fasta of 3 strains
		- \-o <path/to/directory>, --out <path/to/directory>
						Name of output file directory

"""


##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath, existant_file, fasta2dict, directory, parseGFF, loadInList, dict2txt

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
VERSION_DATE='08/07/2016'
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
	parser = argparse.ArgumentParser(prog='buildSNPtoFasta.py', description='''This Programme parse GFF, TAB to add SNP in sequences''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display buildSNPtoFasta.py version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-g', '--gff', metavar="<filename>", type=existant_file, required=True, dest = 'gffFile', help = 'gff files with annotation')
	filesreq.add_argument('-l', '--list', metavar="<filename>", type=existant_file, required=True, dest = 'listKeepFile', help = 'File with geneID to keep')
	filesreq.add_argument('-t', '--tab', metavar="<filename>", type=existant_file, required=True, dest = 'tabFile', help = 'File with SNP')
	filesreq.add_argument('-f', '--fasta', metavar="<path/to/directory>", type = directory, required=True, dest = 'fastaPath', help = 'Directory with fasta of 3 strains')
	filesreq.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'pathOut', help = 'Name of output file directory')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#          Welcome in buildSNPtoFasta (Version " + version + ")          #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les arguments
	gffFile = relativeToAbsolutePath(args.gffFile)
	listKeepFile = relativeToAbsolutePath(args.listKeepFile)
	tabFile = relativeToAbsolutePath(args.tabFile)
	fastaPath = args.fastaPath
	pathFileOut = args.pathOut



	#fastaFile = relativeToAbsolutePath(args.fastaFile)

	print("\t - Input GFF is: %s" % gffFile)
	print("\t - Input listKeppFile is: %s" % listKeepFile)
	print("\t - Input tabFile is: %s" % tabFile)
	print("\t - Input fasta files is: %s" % fastaPath.pathDirectory)
	print("\t - Output fasta files is: %s" % pathFileOut.pathDirectory)

	listKeepID = [ID.replace("Mycfi_gene","gene_") for ID in loadInList(listKeepFile)]
	print("\nTotal ID Keep: %d" % len(listKeepID))


	objGFF = parseGFF(gffFile)
	recordCount = 0
	dicoGenesKeepPosOnScaff = {}
	keepValidList = []
	geneIDsens={}

	for record in objGFF.parseGFF3():

		if record.type == "mRNA" :
			transcriptID = record.attributes["transcriptId"]
			geneID = "gene_"+record.attributes["transcriptId"]


		if record.type == "CDS" and geneID in listKeepID :
			sensDNA = record.strand
			geneIDsens[geneID] = sensDNA

			if geneID not in keepValidList:
				keepValidList.append(geneID)
				recordCount += 1
			for posIter in range(record.start, record.end):
				if record.seqid not in dicoGenesKeepPosOnScaff.keys():
					dicoGenesKeepPosOnScaff[record.seqid] = {posIter:geneID}
				else:
					dicoGenesKeepPosOnScaff[record.seqid].update({posIter:geneID})

	print("\nTotal records: %d" % recordCount)

	print(dict2txt(geneIDsens))



	# Parcours du fichier de SNP
	nblignetotal,ctr = 74215263,0
	dicoSeqBuild = {}
	with open(tabFile, "r") as tabFileRead:
		header = tabFileRead.readline().rstrip().split("\t")
		soucheIndice = header[3:]
		#print(header)

		for line in tabFileRead:

			if ((ctr % 10000 == 0) and (ctr != 0)) or (float(ctr) == nblignetotal):
				percent = (float(ctr)/float(nblignetotal))*100
				sys.stdout.write("\rProcessed up to %0.2f %%..." % percent)
				sys.stdout.flush()

			tabLine = line.rstrip().split("\t")
			scaff = tabLine[0]
			position = int(tabLine[1])
			ref = tabLine[2]
			if scaff in dicoGenesKeepPosOnScaff.keys() and position in dicoGenesKeepPosOnScaff[scaff].keys():
				geneID = dicoGenesKeepPosOnScaff[scaff][position]
				if geneID not in dicoSeqBuild.keys():
					dicoSeqBuild[geneID] = {}


				correctSNP = [snp.replace("R",ref).replace("U","N").replace("F","N") for snp in tabLine[3:]]
				i=0
				for SNP in correctSNP:
					souche = soucheIndice[i]
					#print(souche)
					if souche not in dicoSeqBuild[geneID].keys():
						dicoSeqBuild[geneID][souche] = SNP
					else:
						dicoSeqBuild[geneID][souche] += SNP
					i+=1

			ctr+=1


	#print(dict2txt(dicoSeqBuild))

	# chargement des sequences 3 souches
	for files in fastaPath.lsExtInDirToList("fasta"):
		geneID = files.split("_")[1].replace("gene","gene_")
		dico = fasta2dict(files)
		if geneID in dicoSeqBuild.keys():
			dicoSeqBuild[geneID].update(dico)


	countOnlyN=0
	listSeqNFind = []
	for geneID, dico in dicoSeqBuild.items():

		with open("/work/carlier.j/globalPopGenomicMF/buildSeq62/test/"+geneID+".fasta", "w") as output_handle:
			seqNfind = False
			for souche, txtseq in dico.items():

				if type(txtseq) is not SeqRecord:
					txtseqNdel = txtseq.replace("N","")
					if len(txtseqNdel) == 0:
						countOnlyN += 1
						seqNfind=True
					sensDNA = geneIDsens[geneID]
					if sensDNA == "+":
						record = SeqRecord(Seq(txtseq),id=souche,name=souche, description="MicFi_"+geneID)
					elif sensDNA == "-":
						record = SeqRecord(Seq(txtseq).reverse_complement(),id=souche,name=souche, description="MicFi_"+geneID)
					SeqIO.write(record,output_handle, "fasta")


				else :
					record = SeqRecord(txtseq.seq,id=souche,name=souche, description="MicFi_"+geneID)
					SeqIO.write(record,output_handle, "fasta")


			if seqNfind == True:
				listSeqNFind.append(geneID)


	print("\nNB file with seq Only N: %s" % len(listSeqNFind))
	with open("listToRMgene.txt", "w") as toRMFile:
		txt = "\n".join(listSeqNFind)
		toRMFile.write(txt)



	#print("\n\nExecution summary:")

	#print("  - Outputting \n\
	#Il y a au final %i Sequences gardées sur les %i initial\n\
	#les sequences sont ajoutées dans le fichier %s" %(nbKeep,nbTotal,outputfilename))
	#print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("\n\n#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
