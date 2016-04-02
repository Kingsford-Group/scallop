#!/bin/bash

num=10
dir=./simulation1

for i in `seq 1 $num`
do
	echo "run simulation $i"
	cur=$dir/$i

	./scallop -c ./config -i $cur/expression.gtf -o $cur/refseq.gtf > $cur/refseq.log
#./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop0.gtf -a scallop0 > $cur/scallop0.log
#./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop1.gtf -a scallop1 > $cur/scallop1.log
#./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop2.gtf -a scallop2 > $cur/scallop2.log
#./scallop -c ./config -i $cur/expression.gtf -o $cur/greedy.gtf -a greedy > $cur/greedy.log

done
