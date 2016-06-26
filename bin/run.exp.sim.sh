#!/bin/bash

num=100
dir=./human/exp1

for i in `seq 1 $num`
do
	echo "run simulation $i"
	cur=$dir/$i

	./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop1.gtf -a core -s 5 > $cur/scallop1.log
	./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop2.gtf -a full -s 5 > $cur/scallop2.log
	./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop3.gtf -a greedy -s 5 > $cur/scallop3.log
	./gtfcompare $cur/scallop1.gtf $cur/expression.gtf > $cur/scallop1.cmp
	./gtfcompare $cur/scallop2.gtf $cur/expression.gtf > $cur/scallop2.cmp
	./gtfcompare $cur/scallop3.gtf $cur/expression.gtf > $cur/scallop3.cmp
done
