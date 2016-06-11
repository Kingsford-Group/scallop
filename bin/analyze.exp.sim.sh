#!/bin/bash

num=100
dir=./exp.simulation7

for i in `seq 1 $num`
do
	cur=$dir/$i

	./gtfcompare $cur/scallop1.gtf $cur/expression.gtf > $cur/scallop1.tcmp
	./gtfcompare $cur/scallop2.gtf $cur/expression.gtf > $cur/scallop2.tcmp
	./gtfcompare $cur/scallop3.gtf $cur/expression.gtf > $cur/scallop3.tcmp

	ggg=`cat $cur/scallop2.tcmp | grep summary | cut -f 5 -d " "`
	t2e=`cat $cur/scallop2.tcmp | grep GREATER | wc -l`
	t3e=`cat $cur/scallop3.tcmp | grep GREATER | wc -l`

	t1c=`cat $cur/scallop1.tcmp | grep summary | cut -f 9 -d " "`
	t1d=`cat $cur/scallop1.tcmp | grep summary | cut -f 12 -d " "`

	t2c=`cat $cur/scallop2.tcmp | grep summary | cut -f 9 -d " "`
	t2d=`cat $cur/scallop2.tcmp | grep summary | cut -f 12 -d " "`

	t3c=`cat $cur/scallop3.tcmp | grep summary | cut -f 9 -d " "`
	t3d=`cat $cur/scallop3.tcmp | grep summary | cut -f 12 -d " "`

	echo $i $ggg $t2e $t3e $t1c $t1d $t2c $t2d $t3c $t3d
done
