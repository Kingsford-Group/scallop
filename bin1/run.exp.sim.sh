#!/bin/bash

num=100
dir=./ensembl/$1/exp6

for i in `seq 1 100`
do	
	echo "run simulation $i"

	cur=$dir/$i

	{ time ./scallop -i $cur/expression.gtf -o $cur/scallop1.gtf -a core -s 1 > $cur/scallop1.log; } 2> $cur/scallop1.time
	{ time ./scallop -i $cur/expression.gtf -o $cur/scallop2.gtf -a full -s 1 > $cur/scallop2.log; } 2> $cur/scallop2.time
	{ time ./scallop -i $cur/expression.gtf -o $cur/scallop3.gtf -a greedy -s 1 > $cur/scallop3.log; } 2> $cur/scallop3.time
	./gtfcompare $cur/scallop1.gtf $cur/expression.gtf > $cur/scallop1.cmp
	./gtfcompare $cur/scallop2.gtf $cur/expression.gtf > $cur/scallop2.cmp
	./gtfcompare $cur/scallop3.gtf $cur/expression.gtf > $cur/scallop3.cmp
done

