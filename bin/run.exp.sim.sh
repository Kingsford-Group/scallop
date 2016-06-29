#!/bin/bash

num=100
dir=./ensembl/zebrafish/exp1

list="32 6 27 50 38 16 2 39 44 31 92 93"
for i in `echo $list`
do	
	echo "run simulation $i"

	cur=$dir/$i

	./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop1.gtf -a core -s 1 > $cur/scallop1.log
	./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop2.gtf -a full -s 1 > $cur/scallop2.log
	./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop3.gtf -a greedy -s 1 > $cur/scallop3.log
	./gtfcompare $cur/scallop1.gtf $cur/expression.gtf > $cur/scallop1.cmp
	./gtfcompare $cur/scallop2.gtf $cur/expression.gtf > $cur/scallop2.cmp
	./gtfcompare $cur/scallop3.gtf $cur/expression.gtf > $cur/scallop3.cmp
	./gtfformat $cur/expression.gtf $cur/expression2.cmp
done

