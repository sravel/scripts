#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
## @package struc2runClumpak.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import directory, relativeToAbsolutePath, extant_file, dict2txt

## Python modules
import argparse
from time import localtime, strftime
import re
from subprocess import check_output, STDOUT
##################################################
## Variables Globales
version="0.1"
VERSION_DATE='20/05/2015'

PathPerl5libSeb="/homedir/sravel/perl5/lib/perl5"

drawparams = """
FICHIER MODIF SEB dans Racine
PARAMETERS FOR THE PROGRAM distruct.  YOU WILL NEED TO SET THESE
IN ORDER TO RUN THE PROGRAM.

"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean
        (1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes)
"(d)"   means that this is a double (a real number).

Data settings

#define INFILE_POPQ        casia.popq      // (str) input file of population q's
#define INFILE_INDIVQ      casia.indivq    // (str) input file of individual q's
#define INFILE_LABEL_BELOW casia.names     // (str) input file of labels for below figure
#define INFILE_LABEL_ATOP  labelAtop      // (str) input file of labels for atop figure
#define INFILE_CLUST_PERM  casia.perm     // (str) input file of permutation of clusters to print
#define OUTFILE            casia.ps       //(str) name of output file

#define K	5    // (int) number of clusters
#define NUMPOPS 9    // (int) number of pre-defined populations
#define NUMINDS 210  // (int) number of individuals

Main usage options

#define PRINT_INDIVS      1  // (B) 1 if indiv q's are to be printed, 0 if only population q's
#define PRINT_LABEL_ATOP  1  // (B) print labels above figure
#define PRINT_LABEL_BELOW 1  // (B) print labels below figure
#define PRINT_SEP         0  // (B) print lines to separate populations

Figure appearance

#define FONTHEIGHT 6	// (d) size of font
#define DIST_ABOVE 2	// (d) distance above plot to place text
#define DIST_BELOW -4	// (d) distance below plot to place text
#define BOXHEIGHT  50	// (d) height of the figure
#define INDIVWIDTH 5	// (d) width of an individual

#define ORIENTATION 0	     // (int) 0 for horizontal orientation (default)
			     //       1 for vertical orientation
			     //	      2 for reverse horizontal orientation
                             //       3 for reverse vertical orientation
#define XORIGIN 20		// (d) lower-left x-coordinate of figure
#define YORIGIN 500		// (d) lower-left y-coordinate of figure
#define XSCALE 1		// (d) scale for x direction
#define YSCALE 1		// (d) scale for y direction
#define ANGLE_LABEL_ATOP 90	// (d) angle for labels atop figure (in [0,180])
#define ANGLE_LABEL_BELOW 90    // (d) angle for labels below figure (in [0,180])
#define LINEWIDTH_RIM  3	// (d) width of "pen" for rim of box
#define LINEWIDTH_SEP 0.3	// (d) width of "pen" for separators between pops and for tics
#define LINEWIDTH_IND 0.3	// (d) width of "pen" used for individuals
#define GRAYSCALE 0	        // (B) use grayscale instead of colors
#define ECHO_DATA 0             // (B) print some of the data to the screen
#define REPRINT_DATA 1          // (B) print the data as a comment in the ps file
#define PRINT_INFILE_NAME 0     // (B) print the name of INFILE_POPQ above the figure
                                //     this option is meant for use only with ORIENTATION=0
#define PRINT_COLOR_BREWER 1    // (B) print ColorBrewer settings in the output file
                                //     this option adds 1689 lines and 104656 bytes to the output
                                //     and is required if using ColorBrewer colors
Command line options:

-d drawparams
-K K
-M NUMPOPS
-N NUMINDS
-p input file (population q's)
-i input file (individual q's)
-a input file (labels atop figure)
-b input file (labels below figure)
-c input file (cluster permutation)
-o output file
"""

color = """
1 blue_purple
2 dark_green
3 dark_orange
4 blue
5 light_green
6 light_purple
7 red2
8 yellow
9 color105
10 color109
11 sea_green
12 color122
13 brown
14 pink
15 grey
"""


def buildLabelFile(labelFile, pathDir):
	"""Build label file and labelAtop file"""
	labelFile = relativeToAbsolutePath(labelFile)
	labelAtopOut = open(pathDir+"labelsAtop", "w")
	labelOut = open(pathDir+"label", "w")
	orderlist = []
	with open(labelFile, "r") as labelFileOpen:
		i=1
		for line in labelFileOpen:
			souche, hote = line.rstrip().split("\t")
			hote = hote.replace(" ","_")
			orderlist.append(souche)
			labelAtopOut.write("%i %s\n" %(i, hote))
			labelOut.write("%i\t%s\n" %(i, souche))
			i += 1
	labelAtopOut.close()
	labelOut.close()
	return pathDir+"label", pathDir+"labelsAtop", orderlist

##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog='struc2runClumpak.py', description='''This Programme take directory with output of structure repository and make file to run Clumpak''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display struc2runClumpak version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory', metavar="<path/to/directory>",type = directory, required=True, dest = 'dirPath', help = 'path of result structure')
	filesreq.add_argument('-c', '--clumpak', metavar="<path/to/directory/clumpak>",type = directory, required=True, dest = 'dirPathClumpak', help = 'path of clumpak directory')
	filesreq.add_argument('-l', '--label', metavar="<filename>",type=extant_file, required=True, dest = 'labelFileParam', help = 'File with LABEL, first column name, second top label info')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-dp', '--drawparams', metavar="<filename>",type=extant_file, required=False, dest = 'drawparamsParam', help = 'Check your own drawparams file')
	files.add_argument('-co', '--color', metavar="<filename>",type=extant_file, required=False, dest = 'colorParam', help = 'File with colors (default 15 color max)')


	# Check parameters
	args = parser.parse_args()


	#Welcome message
	print("#################################################################")
	print("#           Welcome in struc2runClumpak (Version " + version + ")           #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les infos passer en argument

	# création de l'objet directory
	workingObjDir = args.dirPath
	clumpakObjDir = args.dirPathClumpak
	labelFileParam = relativeToAbsolutePath(args.labelFileParam)

	# build Drawparam if not add
	if args.drawparamsParam != None:
		drawparamsFile = relativeToAbsolutePath(args.drawparamsParam)
	else:
		drawparamsFile = relativeToAbsolutePath(workingObjDir.pathDirectory+"drawparams")
		with open(drawparamsFile, "w") as drawparamsFileWrite:
			drawparamsFileWrite.write(drawparams)
	# build color if not add
	if args.colorParam != None:
		colorParamFile = relativeToAbsolutePath(args.colorParam)
	else:
		colorParamFile = relativeToAbsolutePath(workingObjDir.pathDirectory+"colorsfile")
		with open(colorParamFile, "w") as colorParamFileWrite:
			colorParamFileWrite.write(color)



	print("Workink Directory: %s" % workingObjDir.pathDirectory)
	#print(workingObjDir.listDir)

	# compte le nombre de répétition et de population
	tabNbPop = []
	tabNbRep = [repo.split("/")[-1].split("_")[-1] for repo in workingObjDir.listDir if "repetition" in repo ]
	for repo in workingObjDir.listDir:
		if "repetition" in repo:
			workingObjDirPop = directory(repo)
			tabNbPop = [repo.split("/")[-1].split("_")[-1] for repo in workingObjDirPop.listDir ]
			break
	if tabNbPop == []:
		print("ERROR struc2runClumpak : path '%s' not contain folder 'repetition_*' with population sub-directory. Make sure take results product by script 'make_structure_dir.py' " % workingObjDir.pathDirectory)
		exit()
	tabNbRepSort = sorted(tabNbRep, key=int)
	tabNbPopSort = sorted(tabNbPop, key=int)
	nbRep ,nbPop = len(tabNbRep), len(tabNbPop)
	print("There are %i repetition of %i population" % (nbRep, nbPop))


	# build label files
	print("Build label files into %s" %workingObjDir.pathDirectory)
	labelFile, labelAtopFile, orderlist = buildLabelFile(labelFileParam, workingObjDir.pathDirectory)
	print(orderlist)


	for rep in tabNbRepSort:																				# boucle sur le nombre de répétition
		if os.path.exists(workingObjDir.pathDirectory+"/repetition_"+str(rep)):
			for pop in tabNbPopSort:
				if os.path.exists(workingObjDir.pathDirectory+"/repetition_"+str(rep)+"/population_"+str(pop)):
					try:
						# remove other file
						os.remove(workingObjDir.pathDirectory+"/repetition_"+str(rep)+"/population_"+str(pop)+"/mainparams")		####### Remove file
						os.remove(workingObjDir.pathDirectory+"/repetition_"+str(rep)+"/population_"+str(pop)+"/extraparams")		####### Remove file
						os.remove(workingObjDir.pathDirectory+"/repetition_"+str(rep)+"/population_"+str(pop)+"/seed.txt")		####### Remove file
					except FileNotFoundError:
						pass

					# open file result to add pop info in order to run clumpak
					filename = os.listdir(workingObjDir.pathDirectory+"/repetition_"+str(rep)+"/population_"+str(pop)+"/")[0]		####### read file name
					print(filename)

					outfile=open(workingObjDir.pathDirectory+"/repetition_"+str(rep)+"/population_"+str(pop)+"/"+filename+'tmp',"w")	# création d'un fichier temps
					toto=0
					dicoFileInfo = {}
					c=1
					inputFile = open(workingObjDir.pathDirectory+"/repetition_"+str(rep)+"/population_"+str(pop)+"/"+filename,"r")
					for ligne in inputFile:
						#print(ligne)
						lignetab = re.sub('\s{2,}', '\t', ligne)
						lligne=lignetab.rstrip().split("\t")
						if "Estimated Allele Frequencies in each cluster" in ligne or ligne.rstrip() == "":
							if toto == 1 or toto == 2:
								toto=0
								compte=1
								#print(dict2txt(dicoFileInfo))
								for souche in orderlist:
									txtOut = str(compte)+" "+" ".join(dicoFileInfo[souche].split(" ")[0:2])+" "+str(compte)+" "+" ".join(dicoFileInfo[souche].split(" ")[3:])
									#print(souche+"\t\t"+txtOut+"\n")
									outfile.write(txtOut+"\n")
									compte+=1

						if toto == 1:
							lignetab = re.sub('\s{1,}', '\t', ligne)
							lligne=lignetab.rstrip().split("\t")
							if c < 100:
								c+=1
								out=lligne[1]+' '+lligne[2]+' '+lligne[3]
								#print(lligne)
								out+=" "+str(lligne[1])+" "
								out+=' '.join(lligne[4:])

							else:
								#print(lligne)
								out=' '.join(lligne[:3])
								out+=" "+str(lligne[0])+" "
								out+=' '.join(lligne[3:])

							souche = out.rstrip().split(" ")[1]
							keep = " ".join(out.rstrip().split(" ")[1:])
							dicoFileInfo[souche] = keep
						if toto == 2:
							lignetab = re.sub('\s{1}', '\t', ligne)
							lligne=lignetab.rstrip().split("\t")
							souche = lligne[1]
							keep = " ".join(lligne[1:])
							dicoFileInfo[souche] = keep



						elif "        Label (%Miss) :  Inferred clusters" not in ligne.rstrip():
							outfile.write(ligne)

						if "        Label (%Miss) :  Inferred clusters" in ligne.rstrip():
							outfile.write("        Label (%Miss) Pop:  Inferred clusters\n")
							toto=1
						if "        Label (%Miss) Pop:  Inferred clusters" in ligne.rstrip():
							toto=2
							c=0

					outfile.close()
					os.remove(workingObjDir.pathDirectory+"/repetition_"+str(rep)+"/population_"+str(pop)+"/"+filename)
					os.rename(workingObjDir.pathDirectory+"/repetition_"+str(rep)+"/population_"+str(pop)+"/"+filename+'tmp',workingObjDir.pathDirectory+"/repetition_"+str(rep)+"/population_"+str(pop)+"/"+filename)


	# Si tous va bien les fichiers sont Ok
	# zip du répertoire
	nameZip = str(workingObjDir.pathDirectory).split("/")[-2]
	zipFile = workingObjDir.pathDirectory+nameZip+".zip"
	if os.path.exists(zipFile):
		print("\nWARNNING struc2runClumpak : path '%s' already contain file %s.zip\nIt will be remove and rebuild\n" %(workingObjDir.pathDirectory, zipFile))
		os.remove(zipFile)
	cmd = "cd "+workingObjDir.pathDirectory+";zip -r "+zipFile+" repetition_*"
	print("Zip directory repetition into %s with file name: %s.zip with command:\n\t%s\n" %(workingObjDir.pathDirectory, zipFile, cmd))
	stream = check_output(cmd, shell=True).decode("utf-8")

	## build label files
	print("Build label files into %s" %workingObjDir.pathDirectory)
	labelFile, labelAtopFile, orderlist = buildLabelFile(labelFileParam, workingObjDir.pathDirectory)

	# copy labelAtop into CLUMPAK/distruct
	cmd = "cp "+labelAtopFile+" "+clumpakObjDir.pathDirectory+"/distruct/labelsAtop"
	print("copy labelAtop into %s with command:\n\t%s\n" %(clumpakObjDir.pathDirectory, cmd))
	stream = check_output(cmd, shell=True).decode("utf-8")

	# make outuput directory of clumpak:
	nameDirectoryOutput = labelFileParam.split("/")[-1].split(".")[0]+"/"
	outputClumpak = workingObjDir.pathDirectory+nameDirectoryOutput
	if os.path.isdir(outputClumpak) == True:
		print("\nWARNNING struc2runClumpak : path '%s' already contain output of CLUMPAK (%s)\nIt will be remove and rebuild\n" %(workingObjDir.pathDirectory, nameDirectoryOutput))
		cmd = "rm -rf "+outputClumpak
		stream = check_output(cmd, shell=True).decode("utf-8")

	cmd = "mkdir "+outputClumpak
	stream = check_output(cmd, shell=True).decode("utf-8")

	# ADD CLUMPAK to perl5lib
	print("LOAD Environment PATH and PERL5LIB:\n")
	oldPERL5LIB = check_output("echo $PERL5LIB", shell=True).decode("utf-8").rstrip()
	os.environ['PERL5LIB'] += ":"+clumpakObjDir.pathDirectory							# add Clumpak to PERL5LIB
	os.environ['PERL5LIB'] += ":"+PathPerl5libSeb							# add module perl for CLUMPAK to PERL5LIB
	oldPATH = check_output("echo $PATH", shell=True).decode("utf-8").rstrip()
	os.environ['PATH'] += ":"+clumpakObjDir.pathDirectory

	# MAKE chmod 755 to Clumpak
	print("Check if exucutable CLUMPAK:\n")
	cmd = "chmod 755 "+clumpakObjDir.pathDirectory+"/mcl/bin/*"
	print(cmd)
	os.system(cmd)
	cmd = "chmod 755 "+clumpakObjDir.pathDirectory+"/CLUMPP/CLUMPP"
	print(cmd)
	os.system(cmd)
	cmd = "chmod 755 "+clumpakObjDir.pathDirectory+"/distruct/distruct1.1"
	print(cmd)
	os.system(cmd)


	# run clumpak
	cmd = "cd "+clumpakObjDir.pathDirectory+"; perl CLUMPAK.pl \
--id "+nameZip+" \
--dir "+outputClumpak+" \
--file "+zipFile+" \
--colors "+colorParamFile+" \
--drawparams "+drawparamsFile+" \
--labels "+labelFile

	print("\nRun CLUMPAK with command:%s\n\nPlease Wait ....\n" % cmd)

	logClumpack = check_output(cmd,stderr=STDOUT, shell=True).decode("utf-8")
	with open(outputClumpak+"Clumpak.log", "w") as log:
		log.write(logClumpack)


	# run BestKByEvanno
	cmd = "cd "+clumpakObjDir.pathDirectory+";\
	perl BestKByEvanno.pl\
	--id "+nameZip+" \
	--d "+outputClumpak+" \
	--f "+zipFile

	print("Run BestKByEvanno with command:%s\n\nPlease Wait ....\n" % cmd)
	logBestKByEvanno = check_output(cmd,stderr=STDOUT, shell=True).decode("utf-8")
	with open(outputClumpak+"BestKByEvanno.log", "w") as log:
		log.write(logClumpack)

	print("RESTORE Environment PATH and PERL5LIB:\n")
	# Rechange PERL5LIB and PATH
	os.environ['PERL5LIB'] = oldPERL5LIB
	os.environ['PATH'] = oldPATH

	print("\n\nExecution summary:")

	print("  - Outputting \n\
\t- Files created in %s." % nameDirectoryOutput)

	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
