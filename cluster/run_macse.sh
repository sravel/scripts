#!/bin/bash -a
# -*- coding: utf-8 -*-
## @package run_macse.sh
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
	printf "#       Run macse on fasta directory Help ( Version $version )       #\n";
	printf "####################################################################\n";
	printf "
 Input:
	directory with fasta's files
 Output:
	directory with NT
	directory with AA
	directory with fasta IN
	directory with sh
	directory with trash

 Exemple Usage: ./run_macse.sh -f ./fasta -m sebastien.ravel@cirad.fr

 Usage: ./run_macse.sh -f {path/to/fasta} -m obiwankenobi@jedi.force
	options:
		-f {path/to/fasta} = path to fasta
		-m {email} = email to add to qsub job end (not mandatory)

		-h = see help\n\n"
	exit 0
}


##################################################
## Parse command line options.
while getopts f:g:m:h: OPT;
	do case $OPT in
		f)	fasta=$OPTARG;;
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

if [ $fasta != "" ] ; then
	#version
	printf "\033[36m ####################################################################\n";
	printf "\033[36m #           Welcome to Run macse directory ( Version $version )            #\n";
	printf "\033[36m ####################################################################\n";

	##################################################
	## Variables Globales

	pathAnalysis=`readlink -m $(dirname $fasta)`"/"
	fastaPath=`readlink -m $fasta`

	NTPath=$pathAnalysis"NT_ALIGN"
	AAPath=$pathAnalysis"AA_ALIGN"
	SHPath=$pathAnalysis"sh"
	trashPath=$pathAnalysis"trash"

	printf "\033[32m \n Working in directory: "$pathAnalysis
	printf "\033[32m \n Fasta were in directory: \n"$fastaPath

	printf "\033[32m \n Output NT_ALIGN were in directory: "$NTPath
	printf "\033[32m \n Output AA_ALIGN were in directory: "$AAPath
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
	if [ -d $NTPath ]; then
		rm -r $NTPath
		mkdir $NTPath
	else
		mkdir $NTPath
	fi
	if [ -d $AAPath ]; then
		rm -r $AAPath
		mkdir $AAPath
	else
		mkdir $AAPath
	fi

	if [ ! -e $pathAnalysis"/submitQsub.sge" ]; then
		touch $pathAnalysis"/submitQsub.sge";
	else
		rm $pathAnalysis"/submitQsub.sge";
		touch $pathAnalysis"/submitQsub.sge";
	fi

	count=1
	for f in $fastaPath/*.fasta ;
	do
		name=$(basename ${f%%.fasta})
		#echo " "$name

		echo "module load system/java/jdk8" > $SHPath"/"$count"_macse.sh"
		echo "java -jar /homedir/sravel/programme/macse_v1.2.jar -prog alignSequences -seq "$f" -ext_gap_ratio 0.00001 -gap_ext 0.00001" >> $SHPath"/"$count"_macse.sh"
		echo "mv "$fastaPath"/"$name"_macse_NT.fasta "$NTPath"/" >> $SHPath"/"$count"_macse.sh"
		echo "mv "$fastaPath"/"$name"_macse_AA.fasta "$AAPath"/" >> $SHPath"/"$count"_macse.sh"

		let count+=1

	done


	echo '#!/bin/bash

#$ -N macse
#$ -cwd
#$ -V
#$ -e '$trashPath'
#$ -o '$trashPath'
#$ -q long.q
#$ -t 1-'$count'
#$ -tc 400
#$ -S /bin/bash

/bin/bash '$pathAnalysis'/sh/${SGE_TASK_ID}_macse.sh

	'>> $pathAnalysis"/submitQsub.sge"


	chmod 755 $pathAnalysis"/submitQsub.sge"


	# Print final infos
	printf "\n\n You want run Macse for "$count" fasta,
 The script are created .sh for all fasta into "$pathAnalysis"sh,\n
 For run all sub-script in qsub, a submitQsub.sge was created, It lunch programm make:\n"

	printf "\033[35m \tqsub "$pathAnalysis"submitQsub.sge "$cmdMail"\n\n"
	# Print end
	printf "\033[36m ####################################################################\n";
	printf "\033[36m #                        End of execution                          #\n";
	printf "\033[36m ####################################################################\n";

# if arguments empty
else
	echo "\033[31m you select fasta = "$fasta
	echo "\033[31m you select mail = "$mail
	printf "\033[31m \n\n You must inform all the required options !!!!!!!!!!!! \n\n"
	help
fi
