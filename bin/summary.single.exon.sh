#!/bin/bash

tag=B439
ver=`stringtie --version`

list="GSM981256.tophat \
	  GSM981244.tophat \
	  GSM984609.tophat\
	  SRR307911.tophat \
	  SRR387662.tophat \
	  GSM981256.star \
	  GSM981244.star \
	  GSM984609.star \
	  SRR307911.star \
	  SRR387662.star"


for i in `echo $list`
do
	x1=`cat $i/scallop.$tag/cuffcmp.all.roc | head -n 1 | cut -f 7  -d " "`
	x2=`cat $i/scallop.$tag/cuffcmp.all.roc | head -n 1 | cut -f 10 -d " "`
	x3=`cat $i/scallop.$tag/cuffcmp.mul.roc | head -n 1 | cut -f 7  -d " "`
	x4=`cat $i/scallop.$tag/cuffcmp.mul.roc | head -n 1 | cut -f 10 -d " "`
	x5=`echo "$x1 - $x3" | bc`
	x6=`echo "$x2 - $x4" | bc`

	y1=`cat $i/stringtie.$ver/cuffcmp.all.roc | head -n 1 | cut -f 7  -d " "`
	y2=`cat $i/stringtie.$ver/cuffcmp.all.roc | head -n 1 | cut -f 10 -d " "`
	y3=`cat $i/stringtie.$ver/cuffcmp.mul.roc | head -n 1 | cut -f 7  -d " "`
	y4=`cat $i/stringtie.$ver/cuffcmp.mul.roc | head -n 1 | cut -f 10 -d " "`
	y5=`echo "$y1 - $y3" | bc`
	y6=`echo "$y2 - $y4" | bc`

	z1=`cat $i/transcomb/cuffcmp.all.roc | head -n 1 | cut -f 7  -d " "`
	z2=`cat $i/transcomb/cuffcmp.all.roc | head -n 1 | cut -f 10 -d " "`
	z3=`cat $i/transcomb/cuffcmp.mul.roc | head -n 1 | cut -f 7  -d " "`
	z4=`cat $i/transcomb/cuffcmp.mul.roc | head -n 1 | cut -f 10 -d " "`
	z5=`echo "$z1 - $z3" | bc`
	z6=`echo "$z2 - $z4" | bc`

	echo $i $x5 $x6 $y5 $y6 $z5 $z6
done
