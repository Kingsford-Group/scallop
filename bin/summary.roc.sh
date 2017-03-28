#!/bin/bash

if [ "$#" != 2 ]; then
	echo "usage ./summary.roc.sh tag abd"
	exit
fi

tag=$1
abd=$2
roc="mul"


scallop="scallop"."$1"."$abd"
stringtie="stringtie.1.3.2d"."$abd"
transcomb="transcomb"."$abd"

list="GSM981256.tophat \
	  GSM981256.star \
	  GSM981256.hisat \
	  GSM981244.tophat \
	  GSM981244.star \
	  GSM981244.hisat \
	  GSM984609.tophat\
	  GSM984609.star \
	  GSM984609.hisat \
	  SRR307911.tophat \
	  SRR307911.star \
	  SRR307911.hisat \
	  SRR387661.tophat \
	  SRR387661.star \
	  SRR387661.hisat"


list="GSM981256.tophat \
	  GSM981244.tophat \
	  GSM984609.tophat\
	  SRR307911.tophat \
	  SRR387661.tophat \
	  GSM981256.star \
	  GSM981244.star \
	  GSM984609.star \
	  SRR307911.star \
	  SRR387661.star \
	  GSM981256.hisat \
	  GSM981244.hisat \
	  GSM984609.hisat \
	  SRR307911.hisat \
	  SRR387661.hisat \
	  SRR545723.tophat \
	  SRR315323.tophat \
	  SRR307903.tophat \
	  SRR315334.tophat \
	  SRR534307.tophat \
	  SRR545723.star \
	  SRR315323.star \
	  SRR307903.star \
	  SRR315334.star \
	  SRR534307.star \
	  SRR545723.hisat \
	  SRR315323.hisat \
	  SRR307903.hisat \
	  SRR315334.hisat \
	  SRR534307.hisat"


for i in `echo $list`
do
	x1=`cat $i/$scallop/gffmul.stats | grep Intron | grep chain | head -n 1 | awk '{print $4}'`
	x2=`cat $i/$scallop/gffmul.stats | grep Intron | grep chain | head -n 1 | awk '{print $6}'`

	y1=`cat $i/$stringtie/gffmul.stats | grep Intron | grep chain | head -n 1 | awk '{print $4}'`
	y2=`cat $i/$stringtie/gffmul.stats | grep Intron | grep chain | head -n 1 | awk '{print $6}'`

	z1=0
	z2=0
	t=`echo $i | cut -f 2 -d "."`
	if [ ! "$t" == "hisat" ]
	then
		z1=`cat $i/$transcomb/gffmul.stats | grep Intron | grep chain | head -n 1 | awk '{print $4}'`
		z2=`cat $i/$transcomb/gffmul.stats | grep Intron | grep chain | head -n 1 | awk '{print $6}'`
	fi

	echo $i $x1 $y1 $z1 $x2 $y2 $z2
done
