#!/bin/bash -a
# -*- coding: utf-8 -*-
## @package run_proteineOrtho.sh
# @author Sebastien Ravel


version=1.0
path=`pwd`

##################################################
## Fonctions
##################################################
# module help
function help
{
	printf "\033[36m####################################################################\n";
	printf "#       Run ProteinOrtho on fastq directory Help ( Version $version )       #\n";
	printf "####################################################################\n";
	printf "
 Input:
	directory with fasta's files
 Output:
	directory with proteinOrtho Results

 Exemple Usage: ./run_proteineOrtho.sh -f ./CDS -t thread -m sebastien.ravel@cirad.fr -s _ALL

 Usage: ./run_proteineOrtho.sh -f {path/to/CDS} -t 10 -m obiwankenobi@jedi.force -s .toto
	options:
		-f {path/to/CDS} = path to fasta with CDS
		-t {int} = number of threads job cluster
		-m {email} = email to add to qsub job end (not mandatory)
		-s {str} = output Suffix Name (not mandatory)

		-h = see help\n\n"
	exit 0
}


##################################################
## Parse command line options.
while getopts f:t:m:h: OPT;
	do case $OPT in
		f)	fasta=$OPTARG;;
		t)	thread=$OPTARG;;
		m)	mail=$OPTARG;;
		s)	suffix=$OPTARG;;
		h)	help;;
		\?)	help;;
	esac
done

if [ $# -eq 0 ]; then
	help
fi


##################################################
## Main code
##################################################

if [ -z ${mail+x} ]; then
	cmdMail=""
else
	cmdMail="-M $mail -m beas"
fi
if [ -z ${suffix+x} ]; then
	suffix=""
fi

if [ $fasta != "" ] && [ $thread != "" ] ; then
	#version
	printf "\033[36m ####################################################################\n";
	printf "\033[36m #        Welcome to Run ProteineOrtho directory ( Version $version )         #\n";
	printf "\033[36m ####################################################################\n";

	##################################################
	## Variables Globales

	pathAnalysis=`readlink -m $(dirname $fasta)`"/"
	fastaPath=`readlink -m $fasta`


	printf "\033[32m \n Working in directory: "$pathAnalysis
	printf "\033[32m \n thread is: "$thread
	printf "\033[32m \n Fasta were in directory: "$fastaPath"\n\n"
	printf "\033[32m \n Suffix is: "$suffix"\n\n"

	if [ -e $pathAnalysis"run_protheineOrtho"$suffix".sh" ]; then
		rm $pathAnalysis"run_protheineOrtho"$suffix".sh"
	fi


	compteur=0
	list=""
	for f in $fastaPath/*.fasta ;
	do
		((compteur++))
		list=$list" "$f

	done


	echo "proteinortho5.pl -cpus=$thread -p=blastn+ -singles -clean -graph -verbose -blastParameters='' -project=phylogenomique"$suffix" $list" >> $pathAnalysis"run_protheineOrtho"$suffix".sh"

	chmod 755 $pathAnalysis"run_protheineOrtho.sh"

	# Print final infos
	printf "\n\n You want run ProteineOrtho for "$compteur" fasta,
 For run use run_protheineOrtho.sh , lunch programm with:\n"

	printf "\033[35m \n\tmodule load compiler/gcc/4.9.2 bioinfo/ncbi-blast/2.2.30 bioinfo/proteinortho/5.11\n"
	printf "\033[35m \tqsub -V -b Y -N ProteinOtho -cwd -q long.q -pe parallel_smp $thread "$pathAnalysis"run_protheineOrtho.sh\n\n"

	# Print end
	printf "\033[36m ####################################################################\n";
	printf "\033[36m #                        End of execution                          #\n";
	printf "\033[36m ####################################################################\n";

# if arguments empty
else
	echo "\033[31m you select fasta = "$fasta
	echo "\033[31m you select thread = "$thread
	echo "\033[31m you select mail = "$mail
	printf "\033[31m \n\n You must inform all the required options !!!!!!!!!!!! \n\n"
	help
fi
