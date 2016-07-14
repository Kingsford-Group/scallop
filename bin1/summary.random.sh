#!/bin/bash

dir=./randomC

for i in `seq 20 20 200`
do
	a=`cat $dir/1000."$i".50/scallop2.cmp | grep GR | wc -l`
	b=`cat $dir/1000."$i".50/scallop3.cmp | grep GR | wc -l`
	ta1=`cat $dir/1000."$i".50/scallop2.time | grep real | cut -f 2 | cut -f 1 -d "m"`
	ta2=`cat $dir/1000."$i".50/scallop2.time | grep real | cut -f 2 | cut -f 2 -d "m" | cut -f 1 -d "."`
	ta=`echo "$ta1 * 60 + $ta2" | bc`
	tb1=`cat $dir/1000."$i".50/scallop3.time | grep real | cut -f 2 | cut -f 1 -d "m"`
	tb2=`cat $dir/1000."$i".50/scallop3.time | grep real | cut -f 2 | cut -f 2 -d "m" | cut -f 1 -d "."`
	tb=`echo "$tb1 * 60 + $tb2" | bc`

	x=`cat $dir/1000."$i".50/scallop2.cmp | head -n 100 | cut -f 2 -d " "`
	y=`echo $x | sed 's/ /+/g'`
	sa=`echo $y | bc`

	x=`cat $dir/1000."$i".50/scallop3.cmp | head -n 100 | cut -f 2 -d " "`
	y=`echo $x | sed 's/ /+/g'`
	sb=`echo $y | bc`

	echo $i $a $b $ta $tb $cc $sa $sb
done
