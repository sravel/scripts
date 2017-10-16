#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package ratioPerSampleMissHet.py
# @author Sajid Ali, Sebastien Ravel

"""
	The ratioPerSampleMissHet script
	===========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 11/04/2017
	:version: 0.1

	Script description
	------------------

	This Programme calculate missing and heterozygosity ratio (diploid data)

	Example
	-------

	>>> ratioPerSampleMissHet.py -vcf test.vcf

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display ratioPerSampleMissHet.py version number and exit

	Input mandatory infos for running:
		- \-vcf <path/to/vcfFile>, --vcfFile <path/to/vcfFile>
						path to vcf File

	Input infos for running with default values:
		- \-o <filename>, --output <filename>
						Output name of tab file (default input.tab)
"""


##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import directory, relativeToAbsolutePath, existant_file, printCol

## Python modules
import argparse
from time import localtime, strftime

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='11/04/2017'

dicoHets = {
			'G/T' : 'K',
			'T/G' : 'K',
			'A/C' : 'M',
			'C/A' : 'M',
			'C/G' : 'S',
			'G/C' : 'S',
			'A/G' : 'R',
			'G/A' : 'R',
			'A/T' : 'W',
			'T/A' : 'W',
			'C/T' : 'Y',
			'T/C' : 'Y'}

listHet = dicoHets.keys()

##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme calculate missing and heterozygosity ratio (diploid data)''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug', action='store_true', dest='debug', help='enter verbose/debug mode')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-vcf', '--vcfFile', metavar="<path/to/vcfFile>",type = existant_file, required=True, dest = 'vcfFile', help = 'path to vcf File')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-o', '--output', metavar="<filename>", required=False, dest = 'outFileName', help = 'Output name of tab file (default input_ratioMissHet.tab)')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#           Welcome in %s (Version %s)            #" %(__file__, version))
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les infos passer en argument
	vcfFile = args.vcfFile
	outFileName = args.outFileName
	if outFileName == None:
		outFileName = vcfFile.replace(".vcf","_ratioMissHet.tab")

	# resume value to user
	print(" - Intput Info:")
	print("\t - VCF files is: %s" % vcfFile)
	print(" - Output Info:")
	print("\t - Output tab created:  %s\n\n" % outFileName)

	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
	#GT:AD:DP:GQ:PL	1:0,51:51:99:1912,0	1:0,18:18:99:727,0	0:46,2:48:99:0,1898	1:0,24:25:99:963,0	1:1,23:24:99:861,0	1:0,10:10:99:421,0	1:0,32:32:99:1348,0	1:0,153:153:99:6386,0

	dicoSampleGenotypes = {}

	# parcours du vcf
	with open(vcfFile, "r") as vcfFileIn, open(outFileName, "w") as outFile:

		headerOut = "Sample\tAllSites\tTotalNoMiss\tRatioMissing\tRatioHet\n"
		outFile.write(headerOut)

		for line in vcfFileIn:
			line = line.rstrip()

			if line.startswith("#"):
				if "#CHROM" in line:
					tabSamples = line.split("\t")[9:]
			else:
				tabLine = line.split("\t")
				#chrom = tabLine[0]
				#pos = tabLine[1]
				#ref = tabLine[3]
				#alleleAlt = tabLine[4]
				tabGenotypes = [infos.split(":")[0] if infos is not "./." else infos for infos in tabLine[9:]]

				#print(tabGenotypes)

				for sample, geno in zip(tabSamples,tabGenotypes):
					if sample not in dicoSampleGenotypes.keys():
						dicoSampleGenotypes[sample] = [geno]
					else:
						dicoSampleGenotypes[sample].append(geno)


		# fin parcours VCF
		for sample, tabGenotypes in dicoSampleGenotypes.items():

			missingCount = float(tabGenotypes.count("./."))
			hetCount = float(sum([allele in listHet for allele in tabGenotypes]))
			allSites = float(len(tabGenotypes))

			totalNoMiss = allSites - missingCount

			ratioMissing = missingCount/allSites

			ratioHet = hetCount/totalNoMiss

			outFile.write("{}\t{:.0f}\t{:.0f}\t{:.4f}\t{:.4f}\n".format(sample, allSites, totalNoMiss, ratioMissing, ratioHet))
			if args.debug: print("{}\t{:.0f}\t{:.0f}\t{:.4f}\t{:.4f}".format(sample, allSites, totalNoMiss, ratioMissing, ratioHet))



	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
