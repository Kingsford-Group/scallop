#!/bin/bash

list="EP.tophat \
	  EP.star \
	  ER.tophat \
	  ER.star \
	  GSM981256.star \
	  GSM981256.hisat \
	  GSM981244.tophat \
	  GSM981244.star \
	  GSM981244.hisat \
	  GSM984609.star \
	  GSM984609.hisat \
	  SRR307911.star \
	  SRR307911.hisat"


for k in `echo $list`
do
	sc=`cat $k/scallop.$1/cmp1.log | tail -n1 | cut -f 6,9,12,15 -d " "`
	st=`cat $k/stringtie/cmp1.log | tail -n1 | cut -f 6,9,12,15 -d " "`
	echo $k $sc $st | awk '{print $1, $2, $6, $3, $7, $4, $8, $5, $9}'
done

list="GSM981256.star \
	  GSM981256.hisat \
	  GSM981244.star \
	  GSM981244.hisat \
	  GSM984609.star \
	  GSM984609.hisat"

xx1="-6"
xx2="-6"
for k in `echo $list`
do
	sc=`cat $k/scallop.$1/cmp1.log | tail -n1 | cut -f 6,9,12,15 -d " "`
	st=`cat $k/stringtie/cmp1.log | tail -n1 | cut -f 6,9,12,15 -d " "`
	x1=`echo $k $sc $st | awk '{print $4 / $8}'`
	x2=`echo $k $sc $st | awk '{print $5 / $9}'`
	xx1="$xx1+$x1"
	xx2="$xx2+$x2"
done

echo "scale = 2; ($xx1) * 100 / 6.0" | bc
echo "scale = 2; ($xx2) * 100 / 6.0" | bc
