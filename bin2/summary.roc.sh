#!/bin/bash

tag=$1
roc="all"
ver=`stringtie --version`

list="EP.tophat \
	  EP.star \
	  ER.tophat \
	  ER.star \
	  GSM981256.tophat \
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
	  SRR387662.tophat \
	  SRR387662.star \
	  SRR387662.hisat"

list="GSM981256.tophat \
	  GSM981256.star \
	  GSM981244.tophat \
	  GSM981244.star \
	  GSM984609.tophat\
	  GSM984609.star \
	  SRR307911.tophat \
	  SRR307911.star \
	  SRR387662.tophat \
	  SRR387662.star"


for i in `echo $list`
do
	x1=`cat $i/scallop.$tag/cuffcmp.$roc.roc | head -n 1 | cut -f 13  -d " "`
	x2=`cat $i/scallop.$tag/cuffcmp.$roc.roc | head -n 1 | cut -f 16 -d " "`

	y1=`cat $i/stringtie.$ver/cuffcmp.$roc.roc | head -n 1 | cut -f 13  -d " "`
	y2=`cat $i/stringtie.$ver/cuffcmp.$roc.roc | head -n 1 | cut -f 16 -d " "`

#echo $i $x1 $y1 $x2 $y2

	z1=0
	z2=0
	t=`echo $i | cut -f 2 -d "."`
	if [ ! "$t" == "hisat" ]
	then
		z1=`cat $i/transcomb/cuffcmp.$roc.roc | head -n 1 | cut -f 13  -d " "`
		z2=`cat $i/transcomb/cuffcmp.$roc.roc | head -n 1 | cut -f 16 -d " "`
	fi

	echo $i $x1 $y1 $z1 $x2 $y2 $z2
done
