#!/bin/bash

version=1.0
runpath=`pwd`
echo $runpath

if [ ! -d $runpath"/sh" ]; then
	mkdir $runpath"/sh";
else
	rm -r $runpath"/sh";
	mkdir $runpath"/sh";
fi
if [ ! -d $runpath"/trash" ]; then
	mkdir $runpath"/trash";
else
	rm -r $runpath"/trash";
	mkdir $runpath"/trash";
fi

if [ ! -e $runpath"/submitQsub.sge" ]; then
	touch $runpath"/submitQsub.sge";
else
	rm $runpath"/submitQsub.sge";
	touch $runpath"/submitQsub.sge";
fi

count=1
for f in $runpath/fasta/*.fasta ;
do
	name=$(basename ${f%%.fasta*})
	echo $name
	echo "module load compiler/gcc/4.9.2 bioinfo/ncbi-blast/2.2.30" > $runpath"/sh/"$count"_blastnr.sh"
	echo "blastx -query "$f" -db /work/BANK/biomaj/nr/nr_2016-05-21/flat/nr -outfmt 5 -out "$runpath"/blastX/"$name".xml" >> $runpath"/sh/"$count"_blastnr.sh"

	let count+=1
done


echo "#!/bin/bash

#$ -N blast
#$ -cwd
#$ -V
#$ -e $runpath/trash/
#$ -o $runpath/trash/
#$ -q long.q
#$ -t 1-$count
#$ -tc 200
#$ -S /bin/bash

/bin/bash '$runpath/sh/'\$SGE_TASK_ID'_blastnr.sh'

">> $runpath"/submitQsub.sge"


chmod 755 $runpath"/submitQsub.sge"
