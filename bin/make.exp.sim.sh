#!/bin/bash

#run these two lines before 
#./purify.gtf.pl $dir/p0.gtf > $dir/p1.gtf
#./fix.gtf.sh $dir/p1.gtf > $dir/p2.gtf

name=p2
gtf=`pwd`/ensembl/human/gtf/$name.gtf
sgtf=`pwd`/ensembl/human/gtf/"$name"_sorted.gtf
dir=`pwd`/ensembl/human/exp1
num=100

mkdir -p $dir

for i in `seq 1 $num`
do
	cur=$dir/$i
	mkdir -p $cur
	ln -sf $gtf $cur/
	ln -sf $sgtf $cur/
	params=$cur/params

	echo "REF_FILE_NAME	$name.gtf" > $params
	echo "NB_MOLECULES	5000000" >> $params
	echo "EXPRESSION_K	-0.1" >> $params
	echo "GEN_DIR			." >> $params
	echo "LOAD_CODING		true" >> $params
	echo "LOAD_NONCODING	false" >> $params
	echo "PRO_FILE_NAME	profile" >> $params
	echo "LIB_FILE_NAME	libfile" >> $params

	flux-simulator -p $params -x
	./merge.exp.pl $gtf $cur/profile > $cur/expression.gtf
done
