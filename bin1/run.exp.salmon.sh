#!/bin/bash

dir=`pwd`/salmon
out=$dir/scallop4
old=$dir/scallop2
exp=$dir/expression.data2
gtf=$dir/gtf

for i in `ls $old`
do	
	if [ -s $old/$i/expression.gtf ]
	then
		echo "run simulation $i"

		cur=$out/$i
		mkdir -p $cur

		ln -sf $old/$i/expression.gtf $cur/expression.gtf

		#./merge.salmon.exp.pl $dir/p3.gtf $exp/$i > $cur/expression.gtf

		{ time ./scallop -i $cur/expression.gtf -o $cur/scallop1.gtf -a core -s 1 > $cur/scallop1.log; } 2> $cur/scallop1.time
		{ time ./scallop -i $cur/expression.gtf -o $cur/scallop2.gtf -a full -s 1 > $cur/scallop2.log; } 2> $cur/scallop2.time
		{ time ./scallop -i $cur/expression.gtf -o $cur/scallop3.gtf -a greedy -s 1 > $cur/scallop3.log; } 2> $cur/scallop3.time
		./gtfcompare $cur/scallop1.gtf $cur/expression.gtf > $cur/scallop1.cmp
		./gtfcompare $cur/scallop2.gtf $cur/expression.gtf > $cur/scallop2.cmp
		./gtfcompare $cur/scallop3.gtf $cur/expression.gtf > $cur/scallop3.cmp
		#./gtfformat $cur/expression.gtf $cur/expression2.gtf
	fi
done
