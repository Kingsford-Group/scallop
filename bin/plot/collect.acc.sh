#!/bin/bash

verA="scallop.B505"
verB="stringtie"
verC="transcomb"

data="GSM981256 \
	  GSM981244 \
	  GSM984609 \
	  SRR307911 \
	  SRR387661 \
	  SRR307903 \
	  SRR315323 \
	  SRR315334 \
	  SRR534307 \
	  SRR545723"

#list="0.01 1 2.5 5 7.5 10 25 50 75 100"
list="0.01"
tmp=./plot

mkdir -p $tmp

for k in `echo $list`
do
		for i in `echo $data`
		do
			x="$j $i"
			a=`cat ../$i.tophat/$verA.$k/gffmul.acc | head -n 1 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.star/$verA.$k/gffmul.acc | head -n 1 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.hisat/$verA.$k/gffmul.acc | head -n 1 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.tophat/$verB.$k/gffmul.acc | head -n 1 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.star/$verB.$k/gffmul.acc | head -n 1 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.hisat/$verB.$k/gffmul.acc | head -n 1 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.tophat/$verC.$k/gffmul.acc | head -n 1 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.star/$verC.$k/gffmul.acc | head -n 1 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			echo $x >> $tmp/acc1.$k

			x="$j $i"
			a=`cat ../$i.tophat/$verA.$k/gffmul.acc | head -n 2 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.star/$verA.$k/gffmul.acc | head -n 2 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.hisat/$verA.$k/gffmul.acc | head -n 2 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.tophat/$verB.$k/gffmul.acc | head -n 2 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.star/$verB.$k/gffmul.acc | head -n 2 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.hisat/$verB.$k/gffmul.acc | head -n 2 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.tophat/$verC.$k/gffmul.acc | head -n 2 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.star/$verC.$k/gffmul.acc | head -n 2 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			echo $x >> $tmp/acc2.$k

			x="$j $i"
			a=`cat ../$i.tophat/$verA.$k/gffmul.acc | head -n 3 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.star/$verA.$k/gffmul.acc | head -n 3 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.hisat/$verA.$k/gffmul.acc | head -n 3 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.tophat/$verB.$k/gffmul.acc | head -n 3 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.star/$verB.$k/gffmul.acc | head -n 3 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.hisat/$verB.$k/gffmul.acc | head -n 3 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.tophat/$verC.$k/gffmul.acc | head -n 3 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			a=`cat ../$i.star/$verC.$k/gffmul.acc | head -n 3 | tail -n 1 | cut -f 13,16 -d " "`; x="$x $a"
			echo $x >> $tmp/acc3.$k
		done
done
