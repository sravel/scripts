#!/bin/bash -a
# -*- coding: utf-8 -*-
## @package runBeast.sh
# @author Sebastien Ravel

version=1.0


##################################################
## Fonctions
##################################################
# module load
function loadModule()
{
	#export OMP_NUM_THREADS=4
	#source /etc/profile.d/modules.sh
	source /homedir/sravel/programme/ScriptsSEB/scripts/loadBeast.sh

	#module load system/java/jre7 system/java/jdk6 compiler/gcc/4.9.2 bioinfo/beagle/4.0 bioinfo/beagle-lib/20150321 bioinfo/BEAST/1.8.1 mpi/openmpi/1.6.5
	printf "Chargement des module OK \n"
}

# module help
function help
{
	printf "\033[36m####################################################################\n";
	printf "#              Run Beast repeats Help ( Version $version )              #\n";
	printf "####################################################################\n";
	printf "
 Input:
	File XML generate by beauti
 Output:
	n directory of repeats into folder with name of xml file

 Exemple Usage: ./runBeast.sh -r 5 -x beauti.xml

 Usage: ./runBeast.sh -r {int} -x {beauti.xml}
	options:
		-r {int} = number of repeats
		-x {beauti.xml} = .xml file build with beauti

		-h = see help\n\n"
	exit 0
}



##################################################
## Parse command line options.
while getopts i:r:x:h: OPT;
	do case $OPT in
		i)	min=$OPTARG;;
		r)	max=$OPTARG;;
		x)	xml=$OPTARG;;
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


if [ -z ${min+x} ]; then
	min=1
fi

if [ $min != "" ] && [ $max != "" ] && [ $xml != "" ]; then
	#version
	printf "\033[36m ####################################################################\n";
	printf "\033[36m #          Welcome to Run Beast repeats ( Version $version )            #\n";
	printf "\033[36m ####################################################################\n";

	##################################################
	## Variables Globales
	#runpath=`pwd`"/"
	runpath=`readlink -m $(dirname $xml)`
	namexml=$(basename $xml .xml)

	pathAnalysis=$runpath$namexml"/"
	printf "\033[32m \n Working in directory: "$runpath"\n\n"
	printf "\033[32m Project name directory: "$pathAnalysis"\n\n"

	loadModule

	# check if min < max
	if [ $min -gt $max ]; then
		printf "\033[31m ERROR "$min" are > than "$max" !!!!!!\n\n";
		exit 0;
	fi

	# check xml file exist
	if [ ! -e $xml ]; then
		printf "\033[31m ERROR "$xml" not exist !!!, please check the path or file, it must be an absolute path from the root \n\n";
		exit 0;
	fi

	# check previous analyse
	if [ -d $pathAnalysis ]; then
		echo -e "\033[33m Directory \""$pathAnalysis"\" already exist !! do you want to remove previous analysis? (y/n):
		y => remove directory and restart analyse;
		n => exit script";
		read answer;
		while [[ $answer != "y" && $answer != "n" ]]; do
			echo -e "\033[31m Please answer y or n !!!!!!!!!";
			read answer;
		done
		if [ $answer == "n" ];then
			printf "\033[31m Program exit !!!!!\n\n"
			printf "\033[36m ####################################################################\n";
			printf "\033[36m #              End of execution of Run Beast repeats               #\n";
			printf "\033[36m ####################################################################\n";
			exit 0;
		fi
		if [ $answer == "y" ];then
			rm -rf $pathAnalysis
			mkdir $pathAnalysis
		fi
	else
		mkdir $pathAnalysis
	fi

	# check SGE PARAMETERS
	qsubTXT="qsub -b Y -V -q long.q -pe parallel_smp 4 -cwd -e trash/ -o trash/";

	echo -e "\033[32m Current qsub command are \""$qsubTXT"\" do you want to keep? (y/n)";
	read answer;
	while [[ $answer != "y" && $answer != "n" ]]; do
		echo -e "\033[31m Please answer y or n !!!!!!!!!";
		read answer;
	done
	if [ $answer == "n" ];then
		echo -e "\033[32m Enter your qsub command";
		read cmd
		while [[ $answer == "" ]]; do
			echo -e "\033[31m Please not keep empty value  !!!!!!!!!";
			read cmd;
		done
		qsubTXT=$cmd
	fi
	if [ $answer == "y" ];then
		qsubTXT="qsub -b Y -V -q long.q -pe parallel_smp 4 -cwd -e trash/ -o trash/";
	fi


	# Run loop of repeats
	for i in $(seq $min $max);
	do

		mkdir $pathAnalysis"repetition"$i
		mkdir $pathAnalysis"repetition"$i"/trash"

		pathXML=$(readlink -m $xml)

		ln -s $pathXML $pathAnalysis"repetition"$i"/"
		#cp  $pathAnalysis"repetition"$i"/"$namexml".xml"

		echo "mpirun -np 4 beast -beagle "$pathAnalysis"repetition"$i"/"$namexml".xml" > $pathAnalysis"repetition"$i"/runBeastRept"$i".sh"

		echo "cd "$pathAnalysis"repetition"$i"/" >> $pathAnalysis"runAllQsub_Beast.sh"
		echo "chmod 755 runBeastRept"$i".sh" >> $pathAnalysis"runAllQsub_Beast.sh"
		echo $qsubTXT" "$pathAnalysis"repetition"$i"/runBeastRept"$i".sh" >> $pathAnalysis"runAllQsub_Beast.sh"


	done


	# Print final infos
	printf "\n\n You want "$max" repeats,
	The script are created all directory of repeats into "$pathAnalysis",\n
	then copy "$namexml".xml into sub-directory and make a runBeastRept.sh for run this repeat\n
	For run all sub-script in qsub, a runAllQsub_Beast.sh was created, It lunch programm make:\n"

	printf "\033[35m \n\texport OMP_NUM_THREADS=4;module load system/java/jre7 system/java/jdk6 compiler/gcc/4.9.2 bioinfo/beagle/4.0 bioinfo/beagle-lib/20150321 bioinfo/BEAST/1.8.1 mpi/openmpi/1.6.5\n"
	printf "\033[35m \tsh "$pathAnalysis"runAllQsub_Beast.sh\n\n"


	# Print end
	printf "\033[36m ####################################################################\n";
	printf "\033[36m #              End of execution of Run Beast repeats               #\n";
	printf "\033[36m ####################################################################\n";

# if arguments empty
else
	echo "you select min = "$min
	echo "you select max = "$max
	echo "you select xml file = "$xml
	printf "\033[31m \n\n You must inform all the required options !!!!!!!!!!!! \n\n"
	help
fi
