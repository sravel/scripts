#!/usr/local/bioinfo/python/2.7.9_build2/bin/python
# -*- coding: utf-8 -*-
# @package vcfToTab.py
# @author Sebastien Ravel

"""
	The compress script
	=========================

	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 28/02/2017
	:version: 0.1

	Script description
	------------------

	This Programme read vcf file (compress or not) and build filter SNP tab

	Example
	-------

	>>> vcfToTab.py -vcf /work/mthierry/89souches_seter_vcf/vcf/ -o ./

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display compress.py version number and exit

	Input mandatory infos for running:
		- \-vcf <path/to/vcfFileDir>, --vcf <path/to/vcfFileDir>
						path to file VCF (compress with gzip or not)

	Input infos for running with default values:
		- \-d <int>, --depth <int>
						Depth threshold (default = 10)
		- \-r <int>, --ratioAD <int>
						RatioAD threshold (default = 0.9)
		- \-o <path/to/outputDir>, --out <path/to/outputDir>
						Name of output directory
		- \-c, --compress
						gzip output file


"""

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import dict2txt, loadInList, dictDict2txt, printCol, relativeToAbsolutePath, existant_file, directory

## Python modules
import argparse
from time import localtime, strftime
import gzip

try:
	import egglib3 as egglib
except Exception as e:
	printCol.red(str(e))
	printCol.red("Must run :\nmodule load compiler/gcc/4.9.3")
	exit()


##################################################
## Variables Globales
version="0.1"
VERSION_DATE='28/02/2017'


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme read vcf file (compress or not) and build filter SNP tab''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-vcf', '--vcf', metavar="<path/to/vcfFileDir>",type = directory, required=True, dest = 'filesDir', help = 'path to file VCF (compress with gzip or not)')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-d', '--depth', metavar="<int>",type=int, default=10, dest = 'thresholdDepth', help = 'Depth threshold (default = 10)')
	files.add_argument('-r', '--ratioAD', metavar="<int>",type=float, default=0.9, dest = 'thresholdRatio', help = 'RatioAD threshold (default = 0.9)')
	files.add_argument('-o', '--out', metavar="<path/to/outputDir>",type = directory, default="./", dest = 'pathOut', help = 'Name of output directory')
	files.add_argument('-c', '--compress',action ='store_true', dest = 'compress', help = 'gzip output file')


	# Check parameters
	args = parser.parse_args()


	#Welcome message
	print("#################################################################")
	print("#              Welcome on vcfToTab (Version " + version + ")               #")
	print("#################################################################\n\n")

	pathFiles = args.filesDir
	pathFilesOut = args.pathOut
	compress = args.compress
	thresholdDepth = args.thresholdDepth
	thresholdRatio = args.thresholdRatio

	# resume value to user
	print(" - Intput Info:")
	print("\t - VCF files are in : %s" % pathFiles.pathDirectory)
	print("\t - %s vcf file count" % len(pathFiles.lsExtInDirToList(["vcf","gz"])))

	print("\t - Depth threshold : %s" % thresholdDepth)
	print("\t - RatioAD threshold : %s" % thresholdRatio)

	print(" - Output Info:")
	if compress:
		print("\t - Output tab will be gzip")
	print("\t - Output tab created into directory:  %s\n\n" % pathFilesOut.pathDirectory)


	for vcfFile in pathFiles.lsExtInDirToList(["vcf","gz"]):
		# récupère le nom du fichier vcf In
		basename = vcfFile.split("/")[-1].split(".")[0]
		extention = vcfFile.split("/")[-1].split(".")[-1]
		outFileName = pathFilesOut.pathDirectory+basename+"_filter_DP"+str(thresholdDepth)+"_ratioAD"+str(thresholdRatio).replace(".","_")+".tab"

		# si le VCF est compresser:
		if extention == "gz":
			import gzip
			vcfFileOpen = gzip.open(vcfFile)
			cache = []
			while True:
				line = vcfFileOpen.readline()
				if line[:2] == '##': cache.append(line)
				elif line[:1] == '#':
					cache.append(line)
					break
				else: raise IOError, 'invalid file'


			header = ''.join(cache)
			# Ouvre le header avec le parseur egglib from_header il sert à comprendre que c'est un fichier vcf
			vcf = egglib.io.VcfParser.from_header(header)
		else:
			# Ouvre le vcf avec le parseur egglib
			vcfFileOpen = egglib.io.VcfParser(str(vcfFile))
			vcf = vcfFileOpen


		samples = [vcf.get_sample(i) for i in xrange(vcf.num_samples)]	# liste des individus

		print 'number of samples: %s\n' % vcf.num_samples
		printCol.green('samples: %s' % "\t".join(samples))


		if compress:
			if ".gz" not in outFileName:
				outputFile = gzip.open(outFileName+'.gz', 'wb')
		else:
			outputFile =  open(outFileName,"w")
		txtout = 'Chromosome\tPosition\tReference\t%s\n' % "\t".join(samples)
		outputFile.write(txtout)


		for line in vcfFileOpen:
			if extention == "gz":							# Parcours des positions du vcf si compresser
				chrom, pos, nall = vcf.read_line(line)
			else: 											# Parcours des positions du vcf si décompresser
				chrom, pos, nall = line

			pos+=1											# egglib commence la position à 0 donc on ajoute 1

			variant = vcf.last_variant()					# récupère la ligne du vcf

			#variant.alleles retourne l'allele de référence en premier puis les autres allèles
			#print "\n\n", nall, variant.alleles
			alt, alt2, alt3 = "", "", ""					# initialise les variants à vide
			if nall == 1:
				ref = variant.alleles[0]
			elif nall == 2:
				ref, alt = variant.alleles					# récupère les allèles dans un tuple
			elif nall == 3:
				ref, alt, alt2 = variant.alleles
			elif nall == 4:
				ref, alt, alt2, alt3 = variant.alleles

			if len(ref) < 2 and len(alt) < 2 and len(alt2) < 2 and len(alt3) < 2:		# test si la présence d'un indel

				try:
					listDepth = []
					listADratio = []
					i = 0
					for dico in variant.samples: 							# variant.samples retourne une liste de dictionnaire / individu (dictionnaire avec 'GT' minimum)
						# On récupère les profondeurs et ratio AD
						if "DP" in dico and dico['DP'] != ():				# dico["DP"] est vide si certains individus n'ont pas de mapping (. dans le vcf)
							#print dico['DP'][0]
							listDepth.append(dico['DP'][0])
						else:
							listDepth.append(0)
						if "AD" in dico and dico["GT"] != ('0',) and dico['AD'] != ():		# test aussi si le génotype n'est pas la reference et si dico["AD"] est vide si certains individus n'ont pas de mapping (. dans le vcf)
							numVariant = dico["GT"][0]										# si plusieurs variants récupére quel est celui de l'individu courent
							nbReadsRef, nbReadsAlt = dico['AD'][0], dico['AD'][int(numVariant)]
							#if pos in [3269,6773, 31671]: printCol.yellow(samples[i]+"  nbReadsRef "+str(nbReadsRef)+"  nbReadsAlt "+str(nbReadsAlt))
							if nbReadsAlt != 0:
								ADratio = nbReadsRef/nbReadsAlt
								listADratio.append(ADratio)
							else:										# si nbReadsAlt = 0 on ajoute 0 pour la valeur du ratio AD
								listADratio.append(0)
						else:											# si pas de Variant pour un individu on ajoute None au ratio AD
							listADratio.append("None")
						i+=1

				except Exception as e:
					print e
					print chrom, pos, samples[i], dico
					exit()


				# Pour récupérer la liste des génotypes de toute nos souches
				genotypes = vcf.get_genotypes()
				#print "genotype",genotypes.as_list()

				# transforme egglib list (0,1) to ref/alt
				variantList = [ str(vartmp[0]).replace("0",ref).replace("1",alt).replace("2",alt2).replace("3",alt3).replace("None","N") for vartmp in genotypes.as_list()[0]]

				try:
					# filtre sur la profondeur et le ratio AD:
					variantFilterList = []
					i=0
					for variant in variantList[0:vcf.num_samples]:
						depth = listDepth[i]
						ratioAD = listADratio[i]
						if int(depth) >= thresholdDepth:
							if ratioAD != "None" and ratioAD <= thresholdRatio :
								variantFilterList.append(variant)
							elif ratioAD != "None" and ratioAD >= thresholdRatio:
								variantFilterList.append("N")
							elif ratioAD == "None" :
								variantFilterList.append(variant)
						else:
							variantFilterList.append("N")
						i+=1
				except Exception as e:
					print e
					print chrom, pos, i, dico
					print listDepth
					print len(listDepth)
					print len(listADratio)
					print len(variantFilterList)
					exit()

				# Ecriture dans le fichier tabuler
				txtout = '%s\t%s\t%s\t%s\n' %(chrom, pos, ref, "\t".join(variantFilterList))
				outputFile.write(txtout)

			#else:
				#print "INDEL at %s %s" % (chrom, pos)
	outputFile.close()

	print("\n\n#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
