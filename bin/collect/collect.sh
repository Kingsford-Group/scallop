#!/bin/bash

dir=B676.1.3.2d
mkdir -p $dir

./match.precision.sh > $dir/match.precision
./match.sensitivity.sh > $dir/match.sensitivity

exit

# collect venn
rm -rf $dir/venn.aligner
rm -rf $dir/venn.algo
tag="W676.1.3.2d"
for i in `cat ../list10`
do
	cat ../$i.$tag/aligner.summary >> $dir/venn.aligner
	cat ../$i.$tag/algo.summary >> $dir/venn.algo
done

for i in `cat ../list10`
do
	echo $i
	./collect.accuracy.sh $i > $dir/$i
done

./collect.quant.sh 1 > $dir/quant1.0.01
./collect.quant.sh 2 > $dir/quant2.0.01
./collect.quant.sh 3 > $dir/quant3.0.01

./collect.class.sh > $dir/class.0.01
./collect.time.sh > $dir/time10
