#!/bin/bash

num=10
dir=./simulation1

for i in `seq 1 $num`
do
	echo "run simulation $i"
	cur=$dir/$i

	./scallop ./config $cur/expression.gtf > $cur/summary
	./scallop ./config $cur/expression.gtf stringtie > $cur/stringtie.log
	./scallop ./config $cur/expression.gtf scallop > $cur/scallop.log

	cat $cur/summary | cut -f 2,3 -d " " | sed 's/,//g' > $cur/reference.path
	cat $cur/stringtie.log | grep solution | cut -f 1,3 -d " " > $cur/stringtie.path
	cat $cur/scallop.log | grep solution | cut -f 1,3 -d " " > $cur/scallop.path
done
