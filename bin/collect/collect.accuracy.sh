#!/bin/bash

if [ "$#" != "1" ]; then
	echo "usage $0 id"
	exit
fi

id=$1

scallop="scallop.B676"
stringtie="stringtie.1.3.2d"
transcomb="transcomb"

list=../list10

# scallop
for aa in `echo "tophat star hisat"`
do
	cc=""
	for abd in `echo "0.01 1 2.5 5 7.5 10 25 50 75 100"`
	do
		x1=`cat ../$id.$aa/$scallop.$abd/gffmul.stats | grep Intron | grep chain | awk '{print $4}'`
		x2=`cat ../$id.$aa/$scallop.$abd/gffmul.stats | grep Intron | grep chain | awk '{print $6}'`
#x1=`cat ../$id.$aa/$scallop.$abd/gffmul.roc | head -n 1 | cut -f 13  -d " "`
#x2=`cat ../$id.$aa/$scallop.$abd/gffmul.roc | head -n 1 | cut -f 16  -d " "`
		cc="$cc$x1 $x2 "
	done
	echo $cc
done

# stringtie
for aa in `echo "tophat star hisat"`
do
	cc=""
	for abd in `echo "0.01 1 2.5 5 7.5 10 25 50 75 100"`
	do
		x1=`cat ../$id.$aa/$stringtie.$abd/gffmul.stats | grep Intron | grep chain | awk '{print $4}'`
		x2=`cat ../$id.$aa/$stringtie.$abd/gffmul.stats | grep Intron | grep chain | awk '{print $6}'`
#x1=`cat ../$id.$aa/$stringtie.$abd/gffmul.roc | head -n 1 | cut -f 13  -d " "`
#x2=`cat ../$id.$aa/$stringtie.$abd/gffmul.roc | head -n 1 | cut -f 16  -d " "`
		cc="$cc$x1 $x2 "
	done
	echo $cc
done

# transcomb
for aa in `echo "tophat star"`
do
	cc=""
	for abd in `echo "0.01 1 2.5 5 7.5 10 25 50 75 100"`
	do
		x1=`cat ../$id.$aa/$transcomb.$abd/gffmul.stats | grep Intron | grep chain | awk '{print $4}'`
		x2=`cat ../$id.$aa/$transcomb.$abd/gffmul.stats | grep Intron | grep chain | awk '{print $6}'`
#x1=`cat ../$id.$aa/$transcomb.$abd/gffmul.roc | head -n 1 | cut -f 13  -d " "`
#x2=`cat ../$id.$aa/$transcomb.$abd/gffmul.roc | head -n 1 | cut -f 16  -d " "`
		cc="$cc$x1 $x2 "
	done
	echo $cc
done
