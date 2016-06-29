#!/bin/bash

num=100
dir=./ensembl/$1/exp1

tmp1=$dir/xx1
tmp2=$dir/xx2

rm -rf $dir/scallop1.cmp
rm -rf $dir/scallop2.cmp
rm -rf $dir/scallop3.cmp
rm -rf $dir/scallop23.cmp
for i in `seq 1 $num`
do
	cur=$dir/$i
	cat $cur/scallop1.cmp | grep -v summary >> $dir/scallop1.cmp
	cat $cur/scallop2.cmp | grep -v summary >> $dir/scallop2.cmp
	cat $cur/scallop3.cmp | grep -v summary >> $dir/scallop3.cmp
#cat $cur/scallop2.cmp | grep -v summary | cut -f 1,2,4 -d " " | sort -k1,1 > $tmp1
#cat $cur/scallop3.cmp | grep -v summary | cut -f 1,2,4 -d " " | sort -k1,1 > $tmp2
#join $tmp1 $tmp2 -11 -21 >> $dir/scallop23.cmp
done

./classify.pl $dir/scallop1.cmp > $dir/summary1
./classify.pl $dir/scallop2.cmp > $dir/summary2
./classify.pl $dir/scallop3.cmp > $dir/summary3
rm -rf $tmp1 $tmp2
