#!/bin/bash

num=10
dir=./simulation1

for i in `seq 1 $num`
do
	cur=$dir/$i

	cat $cur/refseq.log | cut -f 2,3 -d " " | sed 's/,//g' > $cur/refseq.path
	cat $cur/scallop0.log | grep solution | grep assemble0 | cut -f 1,4 -d " " > $cur/scallop0.path
	cat $cur/scallop1.log | grep solution | grep assemble1 | cut -f 1,4 -d " " > $cur/scallop1.path
	cat $cur/scallop2.log | grep solution | grep assemble2 | cut -f 1,4 -d " " > $cur/scallop2.path
	cat $cur/greedy.log   | grep solution | grep greedy    | cut -f 1,4 -d " " > $cur/greedy.path

	./compare.path.pl $cur/scallop0.path $cur/reference.path > $cur/scallop0.pcmp
	./compare.path.pl $cur/scallop1.path $cur/reference.path > $cur/scallop1.pcmp
	./compare.path.pl $cur/scallop2.path $cur/reference.path > $cur/scallop2.pcmp
	./compare.path.pl $cur/greedy.path $cur/reference.path > $cur/greedy.pcmp

	./gtfcompare $cur/scallop0.gtf $cur/refseq.gtf > $cur/scallop0.tcmp
	./gtfcompare $cur/scallop1.gtf $cur/refseq.gtf > $cur/scallop1.tcmp
	./gtfcompare $cur/scallop2.gtf $cur/refseq.gtf > $cur/scallop2.tcmp
	./gtfcompare $cur/greedy.gtf $cur/refseq.gtf > $cur/greedy.tcmp

	total=`cat $cur/refseq.log | wc -l`
	trivial=`cat $cur/refseq.log | grep TRIVIAL | wc -l`
	easy=`cat $cur/refseq.log | grep EASY | wc -l`
	hard=`cat $cur/refseq.log | grep HARD | wc -l`

	ps=`cat $cur/scallop2.pcmp | grep WORSE | wc -l`
	pg=`cat $cur/greedy.pcmp | grep WORSE | wc -l`

	t1a=`cat $cur/scallop1.tcmp | grep summary | cut -f 2 -d " "`
	t1b=`cat $cur/scallop1.tcmp | grep summary | cut -f 5 -d " "`
	t1c=`cat $cur/scallop1.tcmp | grep summary | cut -f 9 -d " "`
	t1d=`cat $cur/scallop1.tcmp | grep summary | cut -f 12 -d " "`

	t2a=`cat $cur/scallop2.tcmp | grep summary | cut -f 2 -d " "`
	t2b=`cat $cur/scallop2.tcmp | grep summary | cut -f 5 -d " "`
	t2c=`cat $cur/scallop2.tcmp | grep summary | cut -f 9 -d " "`
	t2d=`cat $cur/scallop2.tcmp | grep summary | cut -f 12 -d " "`

	ga=`cat $cur/greedy.tcmp | grep summary | cut -f 2 -d " "`
	gb=`cat $cur/greedy.tcmp | grep summary | cut -f 5 -d " "`
	gc=`cat $cur/greedy.tcmp | grep summary | cut -f 9 -d " "`
	gd=`cat $cur/greedy.tcmp | grep summary | cut -f 12 -d " "`

	echo $i $total $trivial $easy $hard $ps $pg $t1a $t1b $t1c $t1d $t2a $t2b $t2c $t2d $ga $gb $gc $gd
done
