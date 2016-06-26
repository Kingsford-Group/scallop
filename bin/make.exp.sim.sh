#!/bin/bash

name=genes1_sorted.gtf
gtf=`pwd`/$name
dir=`pwd`/human/exp1
num=100

mkdir -p $dir

#./purify.gtf.pl $dir/p0.gtf > $dir/p1.gtf
#./fix.gtf.sh $dir/p1.gtf > $dir/p2.gtf

for i in `seq 1 $num`
do
	cur=$dir/$i
	mkdir -p $cur
	ln -sf $gtf $cur/
	params=$cur/params

	echo "REF_FILE_NAME	$name" > $params
	echo "NB_MOLECULES	5000000" >> $params
	echo "EXPRESSION_K	-0.1" >> $params
	echo "GEN_DIR			." >> $params
	echo "LOAD_CODING		true" >> $params
	echo "LOAD_NONCODING	false" >> $params
	echo "PRO_FILE_NAME	profile" >> $params
	echo "LIB_FILE_NAME	libfile" >> $params

	flux-simulator -p $params -x
	
	./merge.exp.pl $gtf $cur/profile > $cur/expression.gtf
	
#	./scallop -c ./config -i $dir/expression.gtf -o $dir/s1.gtf -a core -s 5 > $dir/log1
#	./scallop -c ./config -i $dir/expression.gtf -o $dir/s2.gtf -a full -s 5 > $dir/log2
#	./scallop -c ./config -i $dir/expression.gtf -o $dir/s3.gtf -a greedy -s 5 > $dir/log3
#	
#	./gtfcompare $dir/s1.gtf $dir/expression.gtf > $dir/s1.cmp
#	./gtfcompare $dir/s2.gtf $dir/expression.gtf > $dir/s2.cmp
#	./gtfcompare $dir/s3.gtf $dir/expression.gtf > $dir/s3.cmp
#	
#	./gtfformat $dir/expression.gtf $dir/expression2.gtf
done
