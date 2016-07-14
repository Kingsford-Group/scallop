#!/bin/bash

#dir=./ensembl/$1/exp6
dir=./salmon/scallop4

rm -rf $dir/scallop1.cmp
rm -rf $dir/scallop2.cmp
rm -rf $dir/scallop3.cmp
for i in `ls $dir`
do
	cur=$dir/$i
	cat $cur/scallop1.cmp | grep -v summary >> $dir/scallop1.cmp
	cat $cur/scallop2.cmp | grep -v summary >> $dir/scallop2.cmp
	cat $cur/scallop3.cmp | grep -v summary >> $dir/scallop3.cmp
done

./classify.pl $dir/scallop1.cmp > $dir/summary1
./classify.pl $dir/scallop2.cmp > $dir/summary2
./classify.pl $dir/scallop3.cmp > $dir/summary3
