#!/bin/bash

if [ "$#" != 1 ]; then
	echo "usage ./summary.roc.sh tag"
	exit
fi

tag=$1
roc="mul"

scallop="scallop"."$1"."0.01"
stringtie="stringtie.1.3.2d.0.01"
transcomb="transcomb.0.01"

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


for i in `echo $list`
do
	x1=`cat $i/$scallop/gff$roc.roc | head -n 1 | cut -f 13  -d " "`
	x2=`cat $i/$scallop/gff$roc.roc | head -n 1 | cut -f 16 -d " "`

	y1=`cat $i/$stringtie/gff$roc.roc | head -n 1 | cut -f 13  -d " "`
	y2=`cat $i/$stringtie/gff$roc.roc | head -n 1 | cut -f 16 -d " "`

	z1=0
	z2=0
	t=`echo $i | cut -f 2 -d "."`
	if [ ! "$t" == "hisat" ]
	then
		z1=`cat $i/$transcomb/gff$roc.roc | head -n 1 | cut -f 13  -d " "`
		z2=`cat $i/$transcomb/gff$roc.roc | head -n 1 | cut -f 16 -d " "`
	fi

	echo $i $x1 $y1 $z1 $x2 $y2 $z2
done
