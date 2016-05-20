#!/bin/bash -a
# -*- coding: utf-8 -*-
## @package run_augustus.sh
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
	printf "#       Run augustus on fastq directory Help ( Version $version )       #\n";
	printf "####################################################################\n";
	printf "
 Input:
	directory with fasta's files
 Output:
	directory with GFF

 Exemple Usage: ./run_augustus.sh -f ./fasta -g ./gff -s magnaporthe_grisea -m sebastien.ravel@cirad.fr

 Usage: ./run_augustus.sh -f {path/to/fasta} -g {output/path} -s species -m obiwankenobi@jedi.force
	options:
		-f {path/to/fasta} = path to fasta
		-g {output/path} = path where gff output
		-s {sting} = species name (without space)
		-m {email} = email to add to qsub job end (not mandatory)

		-h = see help\n\n"
	exit 0
}


##################################################
## Parse command line options.
while getopts f:g:s:m:h: OPT;
	do case $OPT in
		f)	fasta=$OPTARG;;
		g)	gff=$OPTARG;;
		s)	species=$OPTARG;;
		m)	mail=$OPTARG;;
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

if [ $fasta != "" ] && [ $gff != "" ] && [ $species != "" ] ; then
	#version
	printf "\033[36m ####################################################################\n";
	printf "\033[36m #        Welcome to Run Augustus directory ( Version $version )          #\n";
	printf "\033[36m ####################################################################\n";

	##################################################
	## Variables Globales

	pathAnalysis=`readlink -m $(dirname $fasta)`"/"
	fastaPath=`readlink -m $fasta`
	gffPath=`readlink -m $gff`

	printf "\033[32m \n Working in directory: "$pathAnalysis"\n"
	printf "\033[32m \n Species is: "$species"\n"
	printf "\033[32m \n Fasta were in directory: "$fastaPath"\n"
	printf "\033[32m \n Output GFF were in directory: "$gffPath"\n\n"

	if [ -d $pathAnalysis"trash" ]; then
		rm -r $pathAnalysis"trash"
		mkdir $pathAnalysis"trash"
	else
		mkdir $pathAnalysis"trash"
	fi
	if [ -d $pathAnalysis"sh" ]; then
		rm -r $pathAnalysis"sh"
		mkdir $pathAnalysis"sh"
	else
		mkdir $pathAnalysis"sh"
	fi
	if [ -e $pathAnalysis"runAllQsub_Augustus.sh" ]; then
		rm $pathAnalysis"runAllQsub_Augustus.sh"
	fi

	compteur=0
	for f in $fastaPath/*.fasta ;
	do
		((compteur++))
		name=$(basename ${f%%.fasta})
		echo " "$name
		echo "augustus --species=$species $f --codingseq=on --outfile=$gffPath/$name.gff" > $pathAnalysis"sh/"$name-augustus.sh
		echo "qsub -N Augustus -b Y -V -q long.q -cwd $cmdMail -e trash/ -o trash/ $pathAnalysis"sh/"$name-augustus.sh" >> $pathAnalysis"runAllQsub_Augustus.sh"

	done


	# Print final infos
	printf "\n\n You want run augustus for "$compteur" fasta,
 The script are created all fasta-augustus.sh for all fasta into "$pathAnalysis"sh,\n
 For run all sub-script in qsub, a runAllQsub_Augustus.sh was created, It lunch programm make:\n"

	printf "\033[35m \n\tmodule load bioinfo/bamtools/8a5d650 compiler/gcc/4.9.2 bioinfo/augustus/3.0.3\n"
	printf "\033[35m \tsh "$pathAnalysis"sh/runAllQsub_Augustus.sh\n\n"


	# Print end
	printf "\033[36m ####################################################################\n";
	printf "\033[36m #              End of execution of Run Beast repeats               #\n";
	printf "\033[36m ####################################################################\n";

# if arguments empty
else
	echo "you select fasta = "$fasta
	echo "you select gff = "$gff
	echo "you select species = "$species
	echo "you select mail = "$mail
	printf "\033[31m \n\n You must inform all the required options !!!!!!!!!!!! \n\n"
	help
fi
