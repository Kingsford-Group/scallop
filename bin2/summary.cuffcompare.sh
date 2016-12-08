#!/bin/bash

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
	  GSM984609.star \
	  GSM984609.hisat \
	  SRR307911.tophat \
	  SRR307911.star \
	  SRR307911.hisat \
	  SRR387662.tophat \
	  SRR387662.star \
	  SRR387662.hisat"

ver=`stringtie --version`

for k in `echo $list`
do
	sc=`./cuffcompare.sh $k/scallop.$1/cuffcmp.stats`
	st=`./cuffcompare.sh $k/stringtie.$ver/cuffcmp.stats`
	tc=`./cuffcompare.sh $k/transcomb/cuffcmp.stats`

	echo $k $sc $st $tc | awk '{print $1, $2, $6, $10, $3, $7, $11, $4, $8, $12, $5, $9, $13}'
done
