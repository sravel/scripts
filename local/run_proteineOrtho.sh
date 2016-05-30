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
	printf "#       Run augustus on fastq directory Help ( Version $version )       #\n";
	printf "####################################################################\n";
	printf "
 Input:
	directory with fasta's files
 Output:
	directory with GFF
	directory with AA
	directory with CDS
	directory with sh
	directory with trash

 Exemple Usage: ./run_proteineOrtho.sh -c ./CDS -t thread -m sebastien.ravel@cirad.fr

 Usage: ./run_proteineOrtho.sh -c {path/to/CDS} -t 10 -m obiwankenobi@jedi.force
	options:
		-c {path/to/CDS} = path to fasta with CDS
		-t {int} = number of threads job cluster
		-m {email} = email to add to qsub job end (not mandatory)

		-h = see help\n\n"
	exit 0
}


##################################################
## Parse command line options.
while getopts f:g:s:m:h: OPT;
	do case $OPT in
		c)	fasta=$OPTARG;;
		t)	thread=$OPTARG;;
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

if [ $fasta != "" ] && [ $thread != "" ] ; then
	#version
	printf "\033[36m ####################################################################\n";
	printf "\033[36m #        Welcome to Run ProteineOrtho directory ( Version $version )         #\n";
	printf "\033[36m ####################################################################\n";

	##################################################
	## Variables Globales

	pathAnalysis=`readlink -m $(dirname $fasta)`"/"
	fastaPath=`readlink -m $fasta`
	gffPath=$pathAnalysis"gff"
	AAPath=$pathAnalysis"AA"
	CDSPath=$pathAnalysis"CDS"
	SHPath=$pathAnalysis"sh"
	trashPath=$pathAnalysis"trash"


#proteinortho5.pl -cpus=3 -p=blastn+ -singles -clean -graph -verbose -blastParameters='' -project=phylogenomique_48souches70-15

	printf "\033[32m \n Working in directory: "$pathAnalysis
	printf "\033[32m \n thread is: "$thread
	printf "\033[32m \n Fasta were in directory: "$fastaPath
	printf "\033[32m \n Output GFF were in directory: "$gffPath
	printf "\033[32m \n Output AA were in directory: "$AAPath
	printf "\033[32m \n Output CDS were in directory: "$CDSPath
	printf "\033[32m \n Output sh were in directory: "$SHPath
	printf "\033[32m \n Output trash were in directory: "$trashPath"\n\n"

	if [ -d $trashPath ]; then
		rm -r $trashPath
		mkdir $trashPath
	else
		mkdir $trashPath
	fi
	if [ -d $SHPath ]; then
		rm -r $SHPath
		mkdir $SHPath
	else
		mkdir $SHPath
	fi
	if [ -d $gffPath ]; then
		rm -r $gffPath
		mkdir $gffPath
	else
		mkdir $gffPath
	fi
	if [ -d $AAPath ]; then
		rm -r $AAPath
		mkdir $AAPath
	else
		mkdir $AAPath
	fi
	if [ -d $CDSPath ]; then
		rm -r $CDSPath
		mkdir $CDSPath
	else
		mkdir $CDSPath
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
		echo "augustus --thread=$thread $f --codingseq=on --protein=on --outfile="$gffPath"/"$name".gff" > $SHPath"/"$name-augustus.sh
		echo "getAnnoFasta.pl "$gffPath"/"$name".gff" >> $SHPath"/"$name-augustus.sh
		echo "mv "$gffPath"/"$name".aa "$AAPath"/" >> $SHPath"/"$name-augustus.sh
		echo "mv "$gffPath"/"$name".codingseq "$CDSPath"/" >> $SHPath"/"$name-augustus.sh
		echo "qsub -N Augustus -b Y -V -q long.q -cwd "$cmdMail" -e "$trashPath" -o "$trashPath" "$SHPath"/"$name"-augustus.sh" >> $pathAnalysis"runAllQsub_Augustus.sh"

		chmod 755 $SHPath"/"$name-augustus.sh

	done

	chmod 755 $pathAnalysis"runAllQsub_Augustus.sh"

	# Print final infos
	printf "\n\n You want run augustus for "$compteur" fasta,
 The script are created all fasta-augustus.sh for all fasta into "$pathAnalysis"sh,\n
 For run all sub-script in qsub, a runAllQsub_Augustus.sh was created, It lunch programm make:\n"

	printf "\033[35m \n\tmodule load compiler/gcc/4.9.2 bioinfo/bamtools/8a5d650 bioinfo/augustus/3.0.3\n"
	printf "\033[35m \tsh "$pathAnalysis"runAllQsub_Augustus.sh\n\n"

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
