#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package make_hapmapfile.py
# @author Lea Picard Sebastien Ravel

"""
	The make_hapmapfile script
	==========================
	:author: Sebastien Ravel, Lea Picard
	:contact: sebastien.ravel@cirad.fr
	:date: 08/07/2016
	:version: 0.1

	Script description
	------------------

	This program takes a SNP matrix with no missing data and returns a hapmap file.\n
	Files will be written in same directory as original file

	Example
	-------

	>>> make_hapmapfile.py -f SNP_table.tab

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display make_hapmapfile.py version number and exit

	Input mandatory infos for running:
		- \-f <filename>, --fileIn <filename>
						Table of SNPs without FU

	Input infos for running with default values:
		- \-o <filename>, --out <filename>
						Name of hapmap file (default outfile.hapmap)
		- \-w <int>, --window <int>
						Minimal window by which SNPs have to be separated, in bp (default = 1 - keep everything)
		- \-c, --chrom
						If used, hapmap files will be produced by chromosomes

"""

##################################################
## Modules
##################################################
## Python modules
import re, argparse, sys, os
sys.path.insert(1,'/homedir/sravel/programme/ScriptsSEB/scripts/modules/')
from MODULES_SEB import relativeToAbsolutePath, existant_file

## Variables Globales
version="0.1"
VERSION_DATE='05/2016'

##################################################
## Main code
##################################################
if __name__ == "__main__":


	# Parameters recovery
	parser = argparse.ArgumentParser(
				prog='make_hapmapfile.py',
				description='''This program takes a SNP matrix with no missing data and returns a hapmap file.
							Files will be written in same directory as original file.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version,
						help='display make_hapmapfile.py version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--fileIn', metavar="<filename>", type = existant_file, required = True,  dest = 'fileIn', help = 'Table of SNPs without FU')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-o', '--out', metavar="<filename>", default="outfile.hapmap", required = False, dest = 'fileOut', help = 'Name of hapmap file (default outfile.hapmap)')
	files.add_argument('-w', '--window', metavar="<int>", type = int, default=1, required = False, dest = 'window', help = 'Minimal window by which SNPs have to be separated, in bp (default = 1 - keep everything)')
	files.add_argument('-c', '--chrom', action='store_true', dest = 'chrom', help = 'If used, hapmap files will be produced by chromosomes')

	# Check parameters
	args = parser.parse_args()

	# get arguments
	fileIn = relativeToAbsolutePath(args.fileIn)
	workdir = "/".join(fileIn.split("/")[:-1])+"/"
	fileOut = workdir+args.fileOut
	window = args.window
	chrom = args.chrom

	if window > 1:
		fileOut = fileOut.split(".")[0]+"_window"+str(window)+".hapmap"
	if window == 0:
		window = 1

	print("#### STARTING")
	print("File will be written in same directory as original file: %s" %workdir)
	print("\n#### CREATING FILES")

	headerHapmap = "rs\talleles\tchrom\tpos\tstrand\tassembly\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t"
	dictFU = {}
	dictSNPScaff = {}
	dictFilesOut = {}

	# by window variables
	oldScaff = ""
	oldPos = 0
	nbSNP = 0

	with open(fileIn) as noFUfile, open(fileOut, "w") as outFile:
		header = noFUfile.readline().rstrip().split("\t")
		headerOut = headerHapmap+"\t".join(header[3:])+"\n"
		outFile.write(headerOut)

		# check numbers
		SNPtot = 0
		rsNumber = 0
		removedLine = 0
		nbLineTot = float(os.popen("wc -l "+fileIn).read().strip("\n").split(" ")[0])
		ctr = 0

		for line in noFUfile:
			line = line.rstrip().split("\t")
			curScaff = line[0]
			SNPtot += 1

			# follow script progression on terminal
			ctr += 1
			if ((ctr % 10000 == 0) and (ctr != 0)) or (float(ctr) == nbLineTot):
				percent = (float(ctr)/float(nbLineTot))*100
				sys.stdout.write("\rProcessed up to %0.2f %%...\t" % percent)
				sys.stdout.flush()

			# check if scaffold name contains a number otherwise attribute 999
			if len(re.findall(r'\d+', curScaff)) == 0:
				print("Scaffold '%s' contains no associated number, automatically attributing 999" %curScaff)
				scaffOut = 999
			else:
				scaffOut = re.findall(r'\d+', curScaff)[0]

			if chrom == True:
				outputName = workdir+fileOut.split('/')[-1].split('.')[0]+"_chr"+scaffOut+".hapmap"

				if scaffOut not in dictFilesOut.keys():
					dictFilesOut[scaffOut] = open(outputName, "w")
					dictFilesOut[scaffOut].write(headerOut)
					#print("\n"+dictFilesOut[scaffOut].name.split("/")[-1])
					sys.stdout.write("\n"+dictFilesOut[scaffOut].name.split("/")[-1])
					sys.stdout.flush()

			# get position of SNP, reference allele and alleles for all samples
			curPos = int(line[1])
			refAll = line[2]
			listAll = line[3:]
			listAll = [allele.replace('R', refAll) for allele in listAll]

			# check number of alleles at position
			dictCount = {x:listAll.count(x) for x in listAll}
			if chrom == True: dictSNPScaff[dictFilesOut[scaffOut]] = (dictSNPScaff.get(dictFilesOut[scaffOut], 0))

			if curScaff != oldScaff:
				oldScaff = curScaff
				oldPos = 0

			diffPos = curPos - oldPos

			if diffPos <= window:
				continue
			else:
				# keep SNP only if bi-allelic
				if len(dictCount.keys()) == 2:
					rsNumber += 1

					if chrom == True: dictSNPScaff[dictFilesOut[scaffOut]] = (dictSNPScaff.get(dictFilesOut[scaffOut], 1))+1

					# create allele couple
					alleleCouple = list(dictCount.keys())

					# create beginning of each line (without genotypes)
					lineHapmap = "rs%i\t%s/%s\t%s\t%s\t+\t2\tNA\tNA\tNA\tNA\tNA\t" %(rsNumber, alleleCouple[0], alleleCouple[1], scaffOut, curPos)
					# pseudo-diploidisation
					listPseudoDiploid = [allele+allele for allele in listAll]
					linePseudoDiploid = "\t".join(listPseudoDiploid)
					# complete line
					lineFileOut = lineHapmap+linePseudoDiploid+"\n"

					if chrom == True: dictFilesOut[scaffOut].write(lineFileOut)
					outFile.write(lineFileOut)

				else:
					removedLine += 1
				oldPos = curPos

		sys.stdout.flush()

		if chrom == True:
			print("\n\n#### CLOSING FILES")
			for keys in sorted(dictFilesOut.keys()):
				openFiles = dictFilesOut[keys]
				print("%s\t(%i SNPs)" %(openFiles.name.split("/")[-1], dictSNPScaff[openFiles]))
				openFiles.close()

		print("\n%i SNPs in the original file." %SNPtot)
		print("%i SNPs successfully written in Hapmap format." %rsNumber)
		print("%i removed SNPs." %removedLine)
