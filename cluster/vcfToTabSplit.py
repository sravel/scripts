#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package vcfToTabSplit.py
# @author Sebastien Ravel

"""
	The vcfToTabSplit script
	===========================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 11/04/2017
	:version: 0.1

	Script description
	------------------

	This Programme build SNP tab with only python

	Example
	-------

	>>> vcfToTabSplit.py -vcf

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display vcfToTabSplit.py version number and exit

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


##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme build SNP tab with only python''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug', action='store_true', dest='debug', help='enter verbose/debug mode')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-vcf', '--vcfFile', metavar="<path/to/vcfFile>",type = existant_file, required=True, dest = 'vcfFile', help = 'path to vcf File')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-o', '--output', metavar="<filename>", required=False, dest = 'outFileName', help = 'Output name of tab file (default input.tab)')

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
		outFileName = vcfFile.replace(".vcf",".tab")

	# resume value to user
	print(" - Intput Info:")
	print("\t - VCF files is: %s" % vcfFile)
	print(" - Output Info:")
	print("\t - Output tab created:  %s\n\n" % outFileName)

	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
	#GT:AD:DP:GQ:PL	1:0,51:51:99:1912,0	1:0,18:18:99:727,0	0:46,2:48:99:0,1898	1:0,24:25:99:963,0	1:1,23:24:99:861,0	1:0,10:10:99:421,0	1:0,32:32:99:1348,0	1:0,153:153:99:6386,0

	# parcours du vcf
	with open(vcfFile, "r") as vcfFileIn, open(outFileName, "w") as outFile:

		for line in vcfFileIn:
			line = line.rstrip()

			if line.startswith("#"):
				if "#CHROM" in line:
					tabSouches = line.split("\t")[9:]
					strSouches = "\t".join(tabSouches)
					print(strSouches)
					txtout = "Chrom\tpos\tref\t%s\n" % strSouches
					outFile.write(txtout)
			else:
				tabLine = line.split("\t")
				chrom = tabLine[0]
				pos = tabLine[1]
				ref = tabLine[3]
				alleleAlt = tabLine[4]
				tabGenotypes = tabLine[9:]

				nallAlt = alleleAlt.count(",")+1
				#print(nallAlt)

				alt, alt2, alt3 = "", "", ""					# initialise les variants à vide
				if nallAlt == 1:
					alt = alleleAlt
				elif nallAlt == 2:
					alt, alt2 = alleleAlt.split(",")					# récupère les allèles dans un tuple
				elif nallAlt == 3:
					alt, alt2, alt3 = alleleAlt.split(",")

				#print(ref, alt, alt2, alt3)

				tabVariant = [ genotype[0].replace("0",ref).replace("1",alt).replace("2",alt2).replace("3",alt3).replace(".","N") for genotype in tabGenotypes ]

				strVariant = "\t".join(tabVariant)
				txtout = "%s\t%s\t%s\t%s\n" % (chrom, pos, ref, strVariant)
				outFile.write(txtout)


	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
