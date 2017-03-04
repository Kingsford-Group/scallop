#!/bin/bash

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
	  SRR387662.tophat \
	  SRR387662.star \
	  SRR387662.hisat"

for x in `echo $list`
do
	id=`echo $x | cut -f 1 -d "."`
	align=`echo $x | cut -f 2 -d "."`
	echo $id $align
	nohup ./run.scallop.sh $x $align B606.0.01 "--library_type first --min_transcript_coverage 0.01" &
done
