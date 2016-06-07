#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
## @package make_structure_dir.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import replace_all, relativeToAbsolutePath, extant_file


## Python modules
import argparse
from time import localtime, strftime

##################################################
## Variables Globales
version="0.2"
VERSION_DATE='05/05/2015'

qsubtxt = "qsub -V -q long.q -cwd "

mainparams = """

KEY PARAMETERS FOR THE PROGRAM structure.  YOU WILL NEED TO SET THESE
IN ORDER TO RUN THE PROGRAM.  VARIOUS OPTIONS CAN BE ADJUSTED IN THE
FILE extraparams.

"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean
        (ie insert 1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes!)

Basic Program Parameters

#define MAXPOPS   **POP**      // (int) number of populations assumed
#define BURNIN    100000   // (int) length of burnin period
#define NUMREPS   300000   // (int) number of MCMC reps after burnin

Input/Output files

#define INFILE   **FILEIN**  // (str) name of input data file
#define OUTFILE  **FILEOUT**  //(str) name of output data file

Data file format

#define NUMINDS    **nIndiv**    // (int) number of diploid individuals in data file
#define NUMLOCI    **nLoci**    // (int) number of loci in data file
#define PLOIDY       1    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line


#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   0     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says
                              whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier

#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data
                             before the genotype data start.

#define MARKERNAMES      1 // (B) data file contains row of marker names
#define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs)
                            // and a row to indicate which alleles are recessive
#define MAPDISTANCES     0  // (B) data file contains row of map distances
                            // between loci


Advanced data file options

#define PHASED           0 // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO        0 // (B) the data for each individual contains a line
                                  indicating phase (linkage model)
#define MARKOVPHASE      0 // (B) the phase info follows a Markov model.
#define NOTAMBIGUOUS  -9 // (int) for use in some analyses of polyploid data



Command line options:

-m mainparams
-e extraparams
-s stratparams
-K MAXPOPS
-L NUMLOCI
-N NUMINDS
-i input file
-o output file
-D SEED
"""


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
#	start=time.clock()
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='make_structure_dir.py', description='''This Programme make arborescence of rep of programme structure''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display make_structure_dir version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-pm', '--popm', metavar="<int>",type = int, required=True, dest = 'nbpopmParam', help = 'Number of pop Max')
	filesreq.add_argument('-i', '--infile', metavar="<fileName>",type=extant_file, required=True, dest = 'inputFile', help = 'input file matrice')
	filesreq.add_argument('-nbi', '--nbIndiv', metavar="<int>",type = int, required=True, dest = 'nbIndivParam', help = 'Number of individus un matrice')
	filesreq.add_argument('-nbm', '--nbMarker', metavar="<int>",type = int, required=True, dest = 'nbMarkerParam', help = 'Number of markers un matrice')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-ri', '--repi', metavar="<int>",type = int, default=1, required=False, dest = 'nbRepiParam', help = 'Number of repetition min (default = 1)')
	files.add_argument('-rm', '--repm', metavar="<int>",type = int, default=10, required=False, dest = 'nbRepmParam', help = 'Number of repetition max (default = 10)')
	files.add_argument('-pi', '--popi', metavar="<int>",type = int, default=1, required=False, dest = 'nbpopiParam', help = 'Number of pop Min (default = 1)')
	files.add_argument('-o', '--outfile', metavar="<PrefixFileName>", required=False, dest = 'outputFile', help = 'output file Prefix (default = name of matrice file)')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#         Welcome in make_structure_dir (Version " + version + ")           #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	nbRepiParam = int(args.nbRepiParam)
	nbRepmParam = int(args.nbRepmParam)
	nbpopiParam = int(args.nbpopiParam)
	nbpopmParam = int(args.nbpopmParam)
	nbIndivParam = args.nbIndivParam
	nbMarkerParam = args.nbMarkerParam
	inputFile = relativeToAbsolutePath(args.inputFile)
	outputFile = args.outputFile



	if outputFile == None:
		outputFile = inputFile.split("/")[-1].split(".")[0]

	workingDir = "/".join(inputFile.split("/")[:-1])+"/"+outputFile+"/"

	print("\t - Input matrice is: %s" % inputFile)
	print("\t - Output prefix name is: %s" % outputFile)
	print("\t - You want %s < K < %s and %s < Repetition < %s" % (nbpopiParam, nbpopmParam, nbRepiParam, nbRepmParam))
	print("\t - Working directory is: %s" % workingDir)

	# ajoute à la variable current_dir le chemin ou est executer le script
	current_dir = os.path.dirname(os.path.abspath(__file__))

	## Test si min < Max
	#if nbRepiParam >= nbRepmParam:
		#print("ERROR: nbRepiParam > nbRepmParam")
		#exit()
	## Test si min < Max
	#if nbpopiParam >= nbpopmParam:
		#print("ERROR: nbpopiParam > nbpopmParam")
		#exit()


	# Test si les sous répertoires existent déjà
	exist=0
	for rep in range(nbRepiParam,nbRepmParam+1):																				# boucle sur le nombre de répétition
		if os.path.exists(workingDir+"/repetition_"+str(rep)):
			print("Warning , folder "+workingDir+"/repetition_"+str(rep)+" already exist !!!!" )
			exist=1
		if exist == 1:
			print("Do you want to remove all analysis? (y/n)\n("+workingDir+"/sh_scripts and /trash will be remove if yes)")
			inp = None
			while inp == None and inp != "y" and inp != "n" and inp != "yes" and inp != "no":
				inp = input()
				if inp == "y":
					os.popen('rm -r '+workingDir+"/repetition_"+str(rep) )
					exist=0
					if os.path.exists(workingDir+"/sh_scripts"):
						os.popen('rm -r '+workingDir+"/sh_scripts" )
					if os.path.exists(workingDir+"/trash"):
						os.popen('rm -r '+workingDir+"/trash" )


	# création des répertoires et fichier mainparams
	os.makedirs(workingDir+"/sh_scripts", exist_ok=True)															# création d'un dossier sh_scripts pour lancer les analyses structures
	os.makedirs(workingDir+"/trash", exist_ok=True)

	shqsub=open(workingDir+"sh_scripts/Qsub_all_structure.sh","w")												# création d'un script qsub pour lancer les sous-scripts

	for rep in range(nbRepiParam,nbRepmParam+1):																	# boucle sur le nombre de répétition
		os.makedirs(workingDir+"/repetition_"+str(rep), exist_ok=True)												# Création du répertoire correspondant
		for pop in range(nbpopiParam,nbpopmParam+1):																# boucle sur la variation du K (np pop)
			#print(str(pop))
			os.makedirs(workingDir+"/repetition_"+str(rep)+"/population_"+str(pop), exist_ok=True)					# création du répertoire de K
			mainparamsOut=open(workingDir+"/repetition_"+str(rep)+"/population_"+str(pop)+"/mainparams","w")		# ouverture du fichier mainparams correspondant
			extraparamsOut=open(workingDir+"/repetition_"+str(rep)+"/population_"+str(pop)+"/extraparams","w")		# ouverture du fichier extraparams correspondant (restera vide)

			# modifie mainarams pour adapter au numero de section , nsection et paths
			dictToReplace = {
			"**POP**"		:	str(pop),
			"**nIndiv**"	:	str(nbIndivParam),
			"**nLoci**"		:	str(nbMarkerParam),
			"**FILEIN**"	:	str(inputFile),
			"**FILEOUT**"	:	outputFile+"_K"+str(pop)+"_run"+str(rep)+".txt"
			}

			mainparamsr = replace_all(dictToReplace, mainparams)

			#mainparamsr = mainparams.replace("**POP**",str(pop)).replace("**FILEIN**",str(inputFile)).replace("**FILEOUT**",outputFile+"_K"+str(pop)+"_run"+str(rep)+".txt")		# modifie mainparams pour adapter au K
			mainparamsOut.write(mainparamsr)																			# Ecrit le mainparams
			mainparamsOut.close()
			extraparamsOut.close()
			#  écriture des scripts pour lancer les annalyses
			shOut=open(workingDir+"/sh_scripts/repetition_"+str(rep)+"_population_"+str(pop)+".sh","w")
			shOut.write("cd "+workingDir+"/repetition_"+str(rep)+"/population_"+str(pop)+"/\n")
			shOut.write("structure")
			shOut.close()
			#  écriture du scripts qsub
			shqsub.write(qsubtxt+workingDir+"/sh_scripts/repetition_"+str(rep)+"_population_"+str(pop)+".sh\n")
	shqsub.close()






		## Display a summary of the execution
	print("\n\nExecution summary:")



	print("  - Outputting \n\
\t- %i directories have been created corresponding to the number repeats, each with %i subdirectories corresponding to the variation number of populations (K).\n\n\
\033[31m  - To launch Structure execute: \n\
\tmodule load bioinfo/structure/2.3.4\n\
\tcd %strash; sh %s\n\
\033[0mon the cluster." %(nbRepmParam, nbpopmParam ,workingDir, shqsub.name))

	print("\033[0m\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
