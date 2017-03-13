#!/bin/bash

dir=B668.1.3.2d
mkdir -p $dir

# collect venn
rm -rf $dir/venn.aligner
rm -rf $dir/venn.algo
tag="W668.1.3.2d"
for i in `cat ../list10`
do
	cat ../$i.$tag/aligner.summary >> $dir/venn.aligner
	cat ../$i.$tag/algo.summary >> $dir/venn.algo
done

exit

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


