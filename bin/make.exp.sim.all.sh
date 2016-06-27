#!/bin/bash

dir=`pwd`/
list="fruitfly chimpanzee mouse rat zebrafish"

for i in `echo $list`
do
	echo $i
	cur=$dir/ensembl/$i/gtf
#./purify.gtf.pl $cur/p0.gtf > $cur/p1.gtf
#./fix.gtf.sh $cur/p1.gtf > $cur/p2.gtf
	flux-simulator -t sortGTF -i $cur/p2.gtf -o $cur/p2_sorted.gtf
	./make.exp.sim.2.sh $i
done
