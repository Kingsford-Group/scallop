#!/bin/bash

tag="B504"

list="GSM981256.star \
	  GSM981256.hisat \
	  GSM981244.star \
	  GSM981244.hisat \
	  GSM984609.star \
	  GSM984609.hisat \
	  SRR307911.star \
	  SRR307911.hisat \
	  SRR387661.star \
	  SRR387661.hisat \
	  GSM981256.tophat \
	  GSM981244.tophat \
	  GSM984609.tophat \
	  SRR307911.tophat \
	  SRR387661.tophat \
	  SRR307903.hisat \
	  SRR307903.star \
	  SRR307903.tophat \
	  SRR315323.hisat \
	  SRR315323.star \
	  SRR315323.tophat \
	  SRR315334.hisat \
	  SRR315334.star \
	  SRR315334.tophat \
	  SRR534307.hisat \
	  SRR534307.star \
	  SRR534307.tophat \
	  SRR545723.hisat \
	  SRR545723.star \
	  SRR545723.tophat"

for i in `echo $list`
do
	id=`echo $i | cut -f 1 -d ","`
	strand=`echo $i |cut -f 2 -d ","`
	aligner=`echo $id | cut -f 2 -d "."`
	a1=`cat $i/scallop.$tag/salmon.quant/logs/salmon_quant.log | grep Mapping | grep rate | awk '{print $8}'`
	a2=`cat $i/stringtie.1.3.1c/salmon.quant/logs/salmon_quant.log | grep Mapping | grep rate | awk '{print $8}'`
	a3=`cat $i/transcomb/salmon.quant/logs/salmon_quant.log | grep Mapping | grep rate | awk '{print $8}'`
	echo $i $a1 $a2 $a3
done
