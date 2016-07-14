#!/bin/bash

dir=`pwd`/
list="worm plant"

for i in `echo $list`
do
	echo $i
	cur=$dir/ensembl/$i/gtf
	./purify.gtf.pl $cur/p0.gtf > $cur/p1.gtf
	./fix.gtf.sh $cur/p1.gtf > $cur/p2.gtf
	./gtfformat $cur/p2.gtf $cur/p3.gtf
	flux-simulator -t sortGTF -i $cur/p2.gtf -o $cur/p2_sorted.gtf
	./make.exp.sim.sh $i
done
