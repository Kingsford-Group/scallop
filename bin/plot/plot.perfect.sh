#!/bin/bash

tag="P3"

dftA="perfect".$tag

tmp=./plot
mkdir -p $tmp

list="GSM981256.ST1 \
	  GSM981244.ST2 \
	  GSM984609.ST3 \
	  SRR387661.TC1 \
	  SRR307911.TC2"

tmp=./plot
mkdir -p $tmp

#echo "REAL.DATA.$measure" | awk '{print toupper($1)}' > $tmp/title
file=$tmp/barplot
rm -rf $file
for ii in `echo $list`
do
	i=`echo $ii | cut -f 1 -d "."`
	j=`echo $ii | cut -f 2 -d "."`
	a0=`cat ../$i.tophat/$dftA/perfect0.roc | head -n 1 | awk '{print $7,$13,$16}'`
	a1=`cat ../$i.tophat/$dftA/perfect1.roc | head -n 1 | awk '{print $7,$13,$16}'`
	a2=`cat ../$i.tophat/$dftA/perfect2.roc | head -n 1 | awk '{print $7,$13,$16}'`
	echo "$i $a0 $a1 $a2" >> $file
done
R CMD BATCH perfect3.R
