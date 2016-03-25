#!/bin/bash

num=10
dir=./simulation1

for i in `seq 1 $num`
do
	cur=$dir/$i

	./compare.path.pl $cur/stringtie.path $cur/reference.path > $cur/stringtie.cmp
	./compare.path.pl $cur/scallop.path $cur/reference.path > $cur/scallop.cmp

	total=`cat $cur/summary | wc -l`
	hard=`cat $cur/summary | grep HARD | wc -l`
	st=`cat $cur/stringtie.cmp | grep WORSE | wc -l`
	sc=`cat $cur/scallop.cmp | grep WORSE | wc -l`

	echo $i $total $hard $st $sc
done
