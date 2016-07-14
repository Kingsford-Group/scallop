#!/bin/bash

num_genes=100
num_exons="1000"
num_transcripts="10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200"
max_length="50"
max_weight=10000
dir=./randomC

for i in `echo $num_exons`
do	
	for j in `echo $num_transcripts`
	do
		for k in `echo $max_length`
		do
			echo "run simulation $i $j $k"
	
			cur=$dir/$i.$j.$k
			mkdir -p $cur
	
			./gtfsimulator $num_genes $i $j $k $max_weight $cur/expression.gtf
	
			{ time ./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop1.gtf -a core -s 1 > $cur/scallop1.log; } 2> $cur/scallop1.time
			{ time ./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop2.gtf -a full -s 1 > $cur/scallop2.log; } 2> $cur/scallop2.time
			{ time ./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop3.gtf -a greedy -s 1 > $cur/scallop3.log; } 2> $cur/scallop3.time
	
			./gtfcompare $cur/scallop1.gtf $cur/expression.gtf > $cur/scallop1.cmp
			./gtfcompare $cur/scallop2.gtf $cur/expression.gtf > $cur/scallop2.cmp
			./gtfcompare $cur/scallop3.gtf $cur/expression.gtf > $cur/scallop3.cmp
		done
	done
done
