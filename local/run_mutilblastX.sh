#!/bin/bash
#!/bin/bash -a
# -*- coding: utf-8 -*-
## @package run_mutilblastX.sh
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
	printf "#       Run MutilblastX on fastq directory Help ( Version $version )       #\n";
	printf "####################################################################\n";
	printf "
 Input:
	directory with fasta's files
 Output:
	directory with blast results
	directory with sh
	directory with trash

 Exemple Usage: ./run_mutilblastX.sh -f ./fasta -o ./blastRes -t blastx -b /work/BANK/biomaj/nr/nr_2016-05-21/flat/nr -m sebastien.ravel@cirad.fr

 Usage: ./run_mutilblastX.sh -f {path/to/fasta} -s {path/to/output} -t blast -b bankBlast -m obiwankenobi@jedi.force
	options:
		*** MANDATORY options:
		-f {path/to/fasta}	= path to fasta
		-o {path/to/output} = path to blast results

		*** Options with default value
		-t {string}			= Blast type (blastx, blastn, ...) (default = blastx)
		-b {path/to/bank}	= path to the databank blast index (default = /work/BANK/biomaj/nr/nr_2016-05-21/flat/nr)
		-f {strint/int}		= outfmt of blast (default = 5)
		-m {email} 			= email to add to qsub job end (not mandatory)

		-h = see help\n\n"
	exit 0
}


##################################################
## Parse command line options.
while getopts f:o:t:m:h: OPT;
	do case $OPT in
		f)	fasta=$OPTARG;;
		o)	out=$OPTARG;;
		t)	type=$OPTARG;;
		b)	bank=$OPTARG;;
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
if [ -z ${bank+x} ]; then
	dbPath="/work/BANK/biomaj/nr/nr_2016-05-21/flat/nr"
else
	dbPath=`readlink -m $bank`
fi
if [ -z ${type+x} ]; then
	typeBlast="blastx"
else
	typeBlast=$type
fi

if [ $fasta != "" ] && [ $out != "" ] ; then
	#version
	printf "\033[36m ####################################################################\n";
	printf "\033[36m #        Welcome to Run MutilblastX directory ( Version $version )         #\n";
	printf "\033[36m ####################################################################\n";

	##################################################
	## Variables Globales


	version=1.0
	runpath=`pwd`
	echo $runpath

	#pathAnalysis=`readlink -m $(dirname $fasta)`"/"
	fastaPath=`readlink -m $fasta`
	outPath=`readlink -m $out`
	pathAnalysis=$outPath"/"

	SHPath=$pathAnalysis"sh"
	trashPath=$pathAnalysis"trash"


	printf "\033[32m \n Working in directory: "$pathAnalysis
	printf "\033[32m \n Fasta were in directory: "$fastaPath
	printf "\033[32m \n blats is: "$typeBlast
	printf "\033[32m \n dataBase is: "$dbPath


	printf "\033[32m \n Output With result Blast were in directory: "$outPath
	printf "\033[32m \n Output sh were in directory: "$SHPath
	printf "\033[32m \n Output trash were in directory: "$trashPath"\n\n"

	if [ -d $SHPath ]; then
		rm -r $SHPath
		mkdir $SHPath
	else
		mkdir $SHPath
	fi
	if [ -d $trashPath ]; then
		rm -r $trashPath
		mkdir $trashPath
	else
		mkdir $trashPath
	fi

	if [ ! -e $pathAnalysis"/submitQsub.sge" ]; then
		touch $pathAnalysis"/submitQsub.sge";
	else
		rm $pathAnalysis"/submitQsub.sge";
		touch $pathAnalysis"/submitQsub.sge";
	fi

	count=1
	for f in $runpath/fasta/*.fasta ;
	do
		name=$(basename ${f%%.fasta*})
		echo $name
		echo "module load compiler/gcc/4.9.2 bioinfo/ncbi-blast/2.2.30" > $pathAnalysis"/sh/"$count"_blastnr.sh"
		echo $typeBlast" -query "$f" -db "$dbPath" -outfmt 5 -out "$pathAnalysis""$name".xml" >> $pathAnalysis"/sh/"$count"_blastnr.sh"

		let count+=1
	done


	echo "#!/bin/bash

	#$ -N blast
	#$ -cwd
	#$ -V
	#$ -e $pathAnalysis/trash/
	#$ -o $pathAnalysis/trash/
	#$ -q long.q
	#$ -t 1-$count
	#$ -tc 200
	#$ -S /bin/bash

	/bin/bash '$pathAnalysis/sh/'\$SGE_TASK_ID'_blastnr.sh'

	">> $pathAnalysis"/submitQsub.sge"


	chmod 755 $pathAnalysis"/submitQsub.sge"


	# Print final infos
	printf "\n\n You want run MutilblastX for "$compteur" fasta,
 The script are created all fasta-MutilblastX.sh for all fasta into "$pathAnalysis"sh,\n
 For run all sub-script in qsub, a submitQsub.sge was created, It lunch programm make:\n"

	printf "\033[35m \tqsub "$pathAnalysis"submitQsub.sge"$cmdMail"\n\n"

	# Print end
	printf "\033[36m ####################################################################\n";
	printf "\033[36m #                        End of execution                          #\n";
	printf "\033[36m ####################################################################\n";

# if arguments empty
else
	echo "\033[31m you select fasta = "$fasta
	echo "\033[31m you select species = "$species
	echo "\033[31m you select mail = "$mail
	printf "\033[31m \n\n You must inform all the required options !!!!!!!!!!!! \n\n"
	help
fi
