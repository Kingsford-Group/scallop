#!/bin/bash

if [ "$#" != "1" ]; then
	echo "usage $0 level(1,2,3)"
	exit
fi

level=$1

scallop="scallop.B676"
stringtie="stringtie.1.3.2d"
transcomb="transcomb"

abd="0.01"
list=../list10

for id in `cat $list`
do
	cc=""

	# scallop
	for aa in `echo "tophat star hisat"`
	do
		x1=`cat ../$id.$aa/$scallop.$abd/gffmul.quant | head -n $level | tail -n 1 | cut -f 13 -d " "`
		x2=`cat ../$id.$aa/$scallop.$abd/gffmul.quant | head -n $level | tail -n 1 | cut -f 16 -d " "`
		cc="$cc$x1 $x2 "
	done

	# stringtie
	for aa in `echo "tophat star hisat"`
	do
		x1=`cat ../$id.$aa/$stringtie.$abd/gffmul.quant | head -n $level | tail -n 1 | cut -f 13 -d " "`
		x2=`cat ../$id.$aa/$stringtie.$abd/gffmul.quant | head -n $level | tail -n 1 | cut -f 16 -d " "`
		cc="$cc$x1 $x2 "
	done

	# transcomb
	for aa in `echo "tophat star"`
	do
		x1=`cat ../$id.$aa/$transcomb.$abd/gffmul.quant | head -n $level | tail -n 1 | cut -f 13 -d " "`
		x2=`cat ../$id.$aa/$transcomb.$abd/gffmul.quant | head -n $level | tail -n 1 | cut -f 16 -d " "`
		cc="$cc$x1 $x2 "
	done

	echo $id $cc
done
