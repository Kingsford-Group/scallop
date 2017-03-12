#!/bin/bash

if [ "$#" != "2" ]
then
	echo "run.sh tag abd"
	exit
fi

tag=$1
abd=$2

#list="SRR534307.hisat"
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

for x in `echo $list`
do
	id=`echo $x | cut -f 1 -d "."`
	align=`echo $x | cut -f 2 -d "."`
	echo $id $align
	nohup ./run.scallop.sh $x $align $tag.$abd "--library_type first --min_transcript_coverage $abd" &
done
