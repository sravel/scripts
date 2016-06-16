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
for f in $runpath/test/*.fas ;
do
	name=$(basename ${f%%.fas*})
	echo $name
	echo "java -jar "$runpath"/macse_v1.01b.jar -prog alignSequences -seq "$f" -ext_gap_ratio 0.00001 -gap_ext 0.00001" > $runpath"/sh/"$count"_alignMacse.sh"


	let count+=1
done



echo "#!/bin/bash

#$ -N MACSE
#$ -cwd
#$ -V
#$ -e $runpath/trash/
#$ -o $runpath/trash/
#$ -q long.q
#$ -t 1-$count
#$ -tc 100
#$ -S /bin/bash

/bin/bash '$runpath/sh/'\$SGE_TASK_ID'_alignMacse.sh'

">> $runpath"/submitQsub.sge"


chmod 755 $runpath"/submitQsub.sge"
