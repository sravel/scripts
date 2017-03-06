#!/usr/local/bioinfo/python/2.7.9_build2/bin/python
# -*- coding: utf-8 -*-
# @package vcfSummaryDepthandRatioAD.py
# @author Sebastien Ravel

"""
	The vcfSummaryDepthandRatioAD script
	====================================

	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 06/03/2017
	:version: 0.1

	Script description
	------------------

	Programme plot histogram of depth and ratio nbReadsRef/nbReadsAlt into one pdf file

	Example
	-------

	>>> vcfSummaryDepthandRatioAD.py -vcf /work/mthierry/89souches_seter_vcf/vcf/souches_publi_teraushi.vcf -o test.pdf

	Help Programm
	-------------

	usage: vcfSummaryDepthandRatioAD.py [-h] [-v] -vcf <path/to/vcfFile> [-o <path/to/pdfFile>]

	optional arguments:
		- \-h, --help
					show this help message and exit
		- \-v, --version
					display vcfSummaryDepthandRatioAD.py version number and exit

	Input mandatory infos for running:
		- \-vcf <path/to/vcfFile>, --vcf <path/to/vcfFile>
					path to file VCF File

	Input infos for running with default values:
		- \-o <path/to/pdfFile>, --out <path/to/pdfFile>
					Name of output pdf file (default = multipage_pdf.pdf)

"""

##################################################
## Modules
##################################################
#Import MODULES_SEB

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import sys, os
from MODULES_SEB import dict2txt, loadInList, dictDict2txt, printCol, relativeToAbsolutePath, existant_file, directory

## Python modules
import argparse
from time import localtime, strftime
import re, subprocess
from subprocess import check_output, STDOUT
try:
	import egglib3 as egglib
except Exception as e:
	printCol.red(str(e))
	printCol.red("Must run :\nmodule load compiler/gcc/4.9.3")
	exit()

## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='24/02/2017'


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme parse vcf file using egglib3 and print depth histogram and ratio nb reads ref/alt for all samples into one pdf ''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-vcf', '--vcf', metavar="<path/to/vcfFile>", required=True, dest = 'vcfFile', help = 'path to file VCF File')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-o', '--out', metavar="<path/to/pdfFile>", default="./multipage_pdf.pdf", dest = 'outFile', help = 'Name of output pdf file (default = multipage_pdf.pdf)')


	# Check parameters
	args = parser.parse_args()


	#Welcome message
	print("#################################################################")
	print("#              Welcome on vcfSummaryDepthandRatioAD (Version " + version + ")               #")
	print("#################################################################\n\n")


	vcfFile = args.vcfFile
	outFile = args.outFile

	# resume value to user
	print(" - Intput Info:")
	print("\t - VCF files is : %s" % vcfFile)

	print(" - Output Info:")
	print("\t - Output pdf is:  %s\n" % outFile)



	# récupère le nom du fichier vcf In
	basename = vcfFile.split("/")[-1].split(".")[0]

	# Ouvre le vcf avec le parseur egglib
	vcf = egglib.io.VcfParser(str(vcfFile))
	samples = [vcf.get_sample(i) for i in xrange(vcf.num_samples)]	# liste des individus

	print 'number of samples: %s\n' % vcf.num_samples
	printCol.green('samples: %s' % "\t".join(samples))

	dicoDepthRatio = {}
	for sample in samples:
		dicoDepthRatio[sample] = {"depth":[], "ratioAD":[]}

	## Parcours des positions du vcf
	nbPrint = 0
	for chrom, pos, nall in vcf:
		pos+=1
		variant = vcf.last_variant()

		#listDepth = []
		#listADratio = []
		i=0
		for dico in variant.samples:							# variant.samples retourne une liste de dictionnaire / individu (dictionnaire avec 'GT' minimum)
			sample = samples[i]
			####### DEPTH
			if "DP" in dico:
				if dico ['DP'] != ():
					dicoDepthRatio[sample]["depth"].append(dico['DP'][0])
					#listDepth.append(dico['DP'][0])
			else:
				#listDepth.append(0)
				dicoDepthRatio[sample]["depth"].append(0)

			####### ratio Ref/Alternate
			if "AD" in dico and dico["GT"] != ('0',) and dico['AD'] != ():		# test aussi si le génotype n'est pas la reference et si dico["AD"] est vide si certains individus n'ont pas de mapping (. dans le vcf)
				numVariant = dico["GT"][0]										# si plusieurs variants récupére quel est celui de l'individu courent
				nbReadsRef, nbReadsAlt = dico['AD'][0], dico['AD'][int(numVariant)]
				if nbReadsAlt != 0:
					ADratio = float(nbReadsRef)/float(nbReadsAlt)
					dicoDepthRatio[sample]["ratioAD"].append(ADratio)
					#listADratio.append(ADratio)
		#i=0
		#for depth, ratioADindiv in zip(listDepth, listADratio):
			#sample = samples[i]
			#dicoDepthRatio[sample]["depth"].append(depth)
			#dicoDepthRatio[sample]["ratioAD"].append(ratioADindiv)
			#i+=1


	# Create the PdfPages object to which we will save the pages:
	with PdfPages(outFile) as pdf:
		for sample, dico in dicoDepthRatio.items():
			print sample
			listDepth = dicoDepthRatio[sample]["depth"]
			listADratio = dicoDepthRatio[sample]["ratioAD"]

			fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=False, figsize=(15,7), facecolor='w', edgecolor='w')
			fig.subplots_adjust(hspace=0.2)

			bins = range(0,260,5)
			ax1.hist(listDepth,bins)										#histogramme des longueurs
			ax1.set_title("%s -- Depth" % sample, size='large', weight='bold')
			ax1.set_xlabel("Coverage", size='large')
			ax1.set_ylabel("Frequency", size='large')
			#plt.close()
			#pdf.savefig()

			bins = np.arange(0,1,0.1)
			ax2.hist(listADratio,bins)										#histogramme des longueurs
			ax2.set_title("%s -- ratio nbReadsRef/nbReadsAlt" % sample, size='large', weight='bold')
			ax2.set_xlabel("nbReadsRef/nbReadsAlt", size='large')
			ax2.set_ylabel("Frequency", size='large')
			pdf.savefig()
			plt.close()


	print("\n\n#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
