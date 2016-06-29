#!/bin/bash

dir=`pwd`/salmon
exp=$dir/expression.data
gtf=$dir/gtf
list=`pwd`/salmon.list1


for i in `cat $list`
do	
	echo "run simulation $i"

	cur=$dir/scallop/$i
	mkdir -p $cur

	./merge.salmon.exp.pl $dir/p3.gtf $exp/$i > $cur/expression.gtf

	./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop1.gtf -a core -s 1 > $cur/scallop1.log
	./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop2.gtf -a full -s 1 > $cur/scallop2.log
	./scallop -c ./config -i $cur/expression.gtf -o $cur/scallop3.gtf -a greedy -s 1 > $cur/scallop3.log
	./gtfcompare $cur/scallop1.gtf $cur/expression.gtf > $cur/scallop1.cmp
	./gtfcompare $cur/scallop2.gtf $cur/expression.gtf > $cur/scallop2.cmp
	./gtfcompare $cur/scallop3.gtf $cur/expression.gtf > $cur/scallop3.cmp
	./gtfformat $cur/expression.gtf $cur/expression2.gtf
done
