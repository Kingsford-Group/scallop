#!/bin/bash

verA="scallop.B505"
verB="stringtie.1.3.2d"
verC="transcomb"
verD="scallop.B511"

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

list="0.01 1 2.5 5 7.5 10 25 50 75 100"

tmp=./results
mkdir -p $tmp

for i in `echo $data`
do
	#verA
	x=""
	for k in `echo $list`
	do
		a=`cat ../$i.tophat/$verA.$k/gffmul.roc | head -n 1 | cut -f 13,16 -d " "`
		x="$x $a"
	done
	echo $x >> $tmp/$i

	x=""
	for k in `echo $list`
	do
		a=`cat ../$i.star/$verA.$k/gffmul.roc | head -n 1 | cut -f 13,16 -d " "`
		x="$x $a"
	done
	echo $x >> $tmp/$i

	x=""
	for k in `echo $list`
	do
		a=`cat ../$i.hisat/$verA.$k/gffmul.roc | head -n 1 | cut -f 13,16 -d " "`
		x="$x $a"
	done
	echo $x >> $tmp/$i


	# verB
	x=""
	for k in `echo $list`
	do
		a=`cat ../$i.tophat/$verB.$k/gffmul.roc | head -n 1 | cut -f 13,16 -d " "`
		x="$x $a"
	done
	echo $x >> $tmp/$i

	x=""
	for k in `echo $list`
	do
		a=`cat ../$i.star/$verB.$k/gffmul.roc | head -n 1 | cut -f 13,16 -d " "`
		x="$x $a"
	done
	echo $x >> $tmp/$i

	x=""
	for k in `echo $list`
	do
		a=`cat ../$i.hisat/$verB.$k/gffmul.roc | head -n 1 | cut -f 13,16 -d " "`
		x="$x $a"
	done
	echo $x >> $tmp/$i

	# verC
	x=""
	for k in `echo $list`
	do
		a=`cat ../$i.tophat/$verC.$k/gffmul.roc | head -n 1 | cut -f 13,16 -d " "`
		x="$x $a"
	done
	echo $x >> $tmp/$i

	x=""
	for k in `echo $list`
	do
		a=`cat ../$i.star/$verC.$k/gffmul.roc | head -n 1 | cut -f 13,16 -d " "`
		x="$x $a"
	done
	echo $x >> $tmp/$i


#	#verD
#	x=""
#	for k in `echo $list`
#	do
#		a=`cat ../$i.tophat/$verD.$k/gffmul.roc | head -n 1 | cut -f 13,16 -d " "`
#		x="$x $a"
#	done
#	echo $x >> $tmp/$i
#
#	x=""
#	for k in `echo $list`
#	do
#		a=`cat ../$i.star/$verD.$k/gffmul.roc | head -n 1 | cut -f 13,16 -d " "`
#		x="$x $a"
#	done
#	echo $x >> $tmp/$i
#
#	x=""
#	for k in `echo $list`
#	do
#		a=`cat ../$i.hisat/$verD.$k/gffmul.roc | head -n 1 | cut -f 13,16 -d " "`
#		x="$x $a"
#	done
#	echo $x >> $tmp/$i

done
