#!/bin/bash

if [ "$#" != "0" ]; then
	echo "usage $0"
	exit
fi

scallop="scallop.B668"
stringtie="stringtie.1.3.2d"
transcomb="transcomb"

abd="0.01"
list=../list10

for exon in `seq 2 20`
do
	for id in `cat $list`
	do
		cc=""

		# scallop
		for aa in `echo "tophat star hisat"`
		do
			x1=`cat ../$id.$aa/$scallop.$abd/gffmul.class | head -n $exon | tail -n 1 | cut -f 6 -d " "`
			x2=`cat ../$id.$aa/$scallop.$abd/gffmul.class | head -n $exon | tail -n 1 | cut -f 9 -d " "`
			cc="$cc$x1 $x2 "
		done

		# stringtie
		for aa in `echo "tophat star hisat"`
		do
			x1=`cat ../$id.$aa/$stringtie.$abd/gffmul.class | head -n $exon | tail -n 1 | cut -f 6 -d " "`
			x2=`cat ../$id.$aa/$stringtie.$abd/gffmul.class | head -n $exon | tail -n 1 | cut -f 9 -d " "`
			cc="$cc$x1 $x2 "
		done

		# transcomb
		for aa in `echo "tophat star"`
		do
			x1=`cat ../$id.$aa/$transcomb.$abd/gffmul.class | head -n $exon | tail -n 1 | cut -f 6 -d " "`
			x2=`cat ../$id.$aa/$transcomb.$abd/gffmul.class | head -n $exon | tail -n 1 | cut -f 9 -d " "`
			cc="$cc$x1 $x2 "
		done

		echo $exon $id $cc
	done
done
