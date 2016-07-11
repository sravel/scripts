#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
## @package run_raxml.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os, shutil
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath, existant_file, directory, printCol

## Python modules
import argparse
from time import localtime, strftime

## BIO Python modules


##################################################
## Variables Globales
version="0.1"
VERSION_DATE='16/06/2016'
debug="False"
#debug="True"


##################################################
## Functions

###################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='run_raxml.py', description='''This Programme run job array to lunch raxml''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display run_raxml.py version number and exit')
	parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesReq = parser.add_argument_group('Input mandatory infos for running')
	filesReq.add_argument('-f', '--fasta', metavar="<path/to/directory/fasta>", type=directory, required=True, dest = 'fastaFileDir', help = 'path to fasta files')
	filesReq.add_argument('-o', '--out', metavar="<path/to/directory>", type = directory, required=True, dest = 'pathOut', help = 'Name of output file directory')

	files = parser.add_argument_group('Input infos for running with default values')
	files.add_argument('-t', '--thread', metavar="<int>",type = int, default=4, required=False, dest = 'nbThreads', help = 'number of threads for raxml (default = 4)')
	files.add_argument('-b', '--bootstrap', metavar="<int>",type = int, default=100, required=False, dest = 'nbBootstrap', help = 'number of nbBootstrap for raxml (default = 100)')
	files.add_argument('-ro', '--raxmloption', metavar="<string>", nargs='*', default=["-f a","-m GTRGAMMA","-x 2","-p 2"],required=False, dest = 'raxmlOptionValue', help = 'Other raxml options (default = "-f a -m GTRGAMMA -x 2 -p 2")')
	files.add_argument('-j', '--nbjob', metavar="<int>", type = int, default=100,required=False, dest = 'nbJobValue', help = 'Number of job array lunch (default = 100)')


	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#              Welcome in run_raxml (Version " + version + ")               #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère les arguments
	pathFastaFile = args.fastaFileDir
	pathFileOut = args.pathOut

	# defaults option
	nbThreads=args.nbThreads
	nbBootstrap=args.nbBootstrap
	raxmlOptionValue=" ".join(args.raxmlOptionValue)

	outputraxmlResDir = pathFileOut.pathDirectory+"raxmlRes/"
	outputSHDir = pathFileOut.pathDirectory+"sh/"
	outputTrashDir = pathFileOut.pathDirectory+"trash/"
	SGENameFile = outputSHDir+"submitQsubraxml.sge"


	# resume value to user
	print(" - Intput Info:")
	print("\t - Working in directory: %s" % pathFileOut.pathDirectory)
	print("\t - Fasta were in directory: %s" % pathFastaFile.pathDirectory)
	print("\t - Number of threads: %s" % nbThreads)
	print("\t - Number of nbBootstrap: %s" % nbBootstrap)
	print("\t - Other options are: %s" % raxmlOptionValue)

	print(" - Output Info:")
	print("\t - Output with result raxml were in directory: %s" % outputraxmlResDir)
	print("\t - Output sh were in directory: %s" % outputSHDir)
	print("\t - Output trash were in directory: %s\n\n" % outputTrashDir)


	# build directory out
	if os.path.exists(outputraxmlResDir):
		printCol.yellow("  Warning , folder %s already exist !!!!" % outputraxmlResDir )
		exist=1
		if exist == 1:
			printCol.yellow("  Do you want to remove all analysis? (y/n)\n  ("+outputSHDir+" and "+outputTrashDir+" will be remove if yes)")
			inp = None
			while inp == None and inp != "y" and inp != "n" and inp != "yes" and inp != "no":
				inp = input()
				if inp == "y":
					#os.popen('rm -rf '+outputraxmlResDir+" && mkdir "+outputraxmlResDir)
					#os.popen('rm -rf '+outputSHDir+" && mkdir "+outputSHDir)
					#os.popen('rm -rf '+outputTrashDir+" && mkdir "+outputTrashDir)
					shutil.rmtree(outputraxmlResDir)
					shutil.rmtree(outputSHDir)
					shutil.rmtree(outputTrashDir)

	if not os.path.exists(outputraxmlResDir):
		os.makedirs(outputraxmlResDir)
		os.makedirs(outputraxmlResDir+"RAxML_bipartitionsBranchLabels")
		os.makedirs(outputraxmlResDir+"RAxML_bipartitions")
		os.makedirs(outputraxmlResDir+"RAxML_bestTree")
		os.makedirs(outputraxmlResDir+"RAxML_bootstrap")
		os.makedirs(outputraxmlResDir+"RAxML_info")
		os.makedirs(outputraxmlResDir+"RAxML_reduced")

	if not os.path.exists(outputSHDir):
		os.makedirs(outputSHDir)								# création d'un dossier sh_scripts pour lancer les analyses structures
	if not os.path.exists(outputTrashDir):
		os.makedirs(outputTrashDir)

	count = 1
	for fasta in pathFastaFile.lsExtInDirToList("fasta"):
		basenameFasta = fasta.split("/")[-1].split(".")[0]

		with open(outputSHDir+str(count)+"_raxml.sh", "w") as shScript:
			shScript.write("module load mpi/openmpi/1.6.5 compiler/gcc/4.9.2 bioinfo/RAxML/8.1.17\n")
			shScript.write('sed -i s/"\!"/"n"/g '+fasta+'\n')
			raxmlcmd = "raxmlHPC-PTHREADS -T %s -#%s -n %s %s -s %s\n" % (nbThreads, nbBootstrap, basenameFasta, raxmlOptionValue, fasta)
			if args.debug == "True" : print(raxmlcmd)
			shScript.write(raxmlcmd)
			shScript.write("mv "+pathFileOut.pathDirectory+"RAxML_bipartitionsBranchLabels."+basenameFasta+" "+outputraxmlResDir+"RAxML_bipartitionsBranchLabels/\n")
			shScript.write("mv "+pathFileOut.pathDirectory+"RAxML_bipartitions."+basenameFasta+" "+outputraxmlResDir+"RAxML_bipartitions/\n")
			shScript.write("mv "+pathFileOut.pathDirectory+"RAxML_bestTree."+basenameFasta+" "+outputraxmlResDir+"RAxML_bestTree/\n")
			shScript.write("mv "+pathFileOut.pathDirectory+"RAxML_bootstrap."+basenameFasta+" "+outputraxmlResDir+"RAxML_bootstrap/\n")
			shScript.write("mv "+pathFileOut.pathDirectory+"RAxML_info."+basenameFasta+" "+outputraxmlResDir+"RAxML_info/\n")
			shScript.write("mv "+fasta+".reduced "+outputraxmlResDir+"RAxML_reduced/\n")




		count+=1

#raxmlHPC-PTHREADS -T 24 -#100 -n 48souchesAllCDS -f a -m GTRGAMMA -x 2 -p 2 -s /work/sravel/align48/raxml/48souchesAllCDS.fasta

#mkdir RAxML_bipartitionsBranchLabels
#mv RAxML_bipartitionsBranchLabels.M* ./RAxML_bipartitionsBranchLabels


#mkdir RAxML_bipartitions
#mv RAxML_bipartitions.M* ./RAxML_bipartitions

#mkdir RAxML_bestTree
#mv RAxML_bestTree.M* ./RAxML_bestTree

#mkdir RAxML_info
#mv RAxML_info.M* ./RAxML_info

	headerSGE = """
#!/bin/bash

#$ -N raxml
#$ -cwd
#$ -V
#$ -e """+outputTrashDir+"""
#$ -o """+outputTrashDir+"""
#$ -q long.q
#$ -pe parallel_smp """+str(nbThreads)+"""
#$ -t 1-"""+str(count-1)+"""
#$ -tc """+str(args.nbJobValue)+"""
#$ -S /bin/bash

/bin/bash """+outputSHDir+"""${SGE_TASK_ID}_raxml.sh"""


	with open(SGENameFile, "w") as SGEFile:
		SGEFile.write(headerSGE)

	print("\n - Execution summary:")

	print("\n  You want run Mutilraxml for %s fasta,\
 The script are created all fasta-Mutilraxml.sh for all fasta into %s,\n\
 For run all sub-script in qsub, %s was created.\n\
 It lunch programm with job array and run %s job max:\n" %(count-1,outputSHDir,SGENameFile, args.nbJobValue))
	printCol.green("\tqsub %s" % SGENameFile)
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
