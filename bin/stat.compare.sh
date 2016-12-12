#!/bin/bash

dir=$1

tmp=./cmp.tmp
rm -rf $tmp

for i in `ls $dir | grep log`
do
	j=${i/.log/}

	b=`cat $dir/$i | grep Bundle | wc -l`

	if [ "$b" -ge 2 ]
	then
		continue
	fi

	cat $dir/$i | grep summary | grep TP >> $tmp
done

a=`cat $tmp | grep splice | cut -f  6 -d " "`; tp1=`echo $a | sed 's/ /+/g' | bc`;
a=`cat $tmp | grep splice | cut -f  9 -d " "`; fp1=`echo $a | sed 's/ /+/g' | bc`;
a=`cat $tmp | grep splice | cut -f 12 -d " "`; fn1=`echo $a | sed 's/ /+/g' | bc`;

a=`cat $tmp | grep boundary | cut -f  6 -d " "`; tp2=`echo $a | sed 's/ /+/g' | bc`;
a=`cat $tmp | grep boundary | cut -f  9 -d " "`; fp2=`echo $a | sed 's/ /+/g' | bc`;
a=`cat $tmp | grep boundary | cut -f 12 -d " "`; fn2=`echo $a | sed 's/ /+/g' | bc`;

echo $1 splice-positions TP = $tp1 FP = $fp1 FN = $fn1 boundaries TP = $tp2 FP = $fp2 FN = $fn2

rm -rf $tmp
