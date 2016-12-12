#!/bin/bash

tag="B466"
measure="all"

dftA="scallop".$tag
dftB="stringtie.1.3.1c"
dftC="transcomb"

tmp=./plot
mkdir -p $tmp

list="GSM981256.ST1 \
	  GSM981244.ST2 \
	  GSM984609.ST3 \
	  SRR387662.TC1 \
	  SRR307911.TC2"

tmp=./plot
mkdir -p $tmp

echo "REAL.DATA.$measure" | awk '{print toupper($1)}' > $tmp/title
file=$tmp/barplot
rm -rf $file
for ii in `echo $list`
do
	i=`echo $ii | cut -f 1 -d "."`
	j=`echo $ii | cut -f 2 -d "."`
	a1=`cat ../$i.tophat/$dftA/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	a2=`cat ../$i.star/$dftA/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	a3=`cat ../$i.hisat/$dftA/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	b1=`cat ../$i.tophat/$dftB/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	b2=`cat ../$i.star/$dftB/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	b3=`cat ../$i.hisat/$dftB/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	c1=`cat ../$i.tophat/$dftC/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	c2=`cat ../$i.star/$dftC/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	c3="0 0 0"
	echo "$j $a1 $a2 $a3 $b1 $b2 $b3 $c1 $c2 $c3" >> $file
done
R CMD BATCH barplot.9.R
mv $tmp/pre.pdf $tmp/sra.$measure.pre.pdf
mv $tmp/sen.pdf $tmp/sra.$measure.sen.pdf
mv $file $tmp/sra.$measure.sen.summary

list="EP ER"

tmp=./plot
mkdir -p $tmp

echo "REAL.DATA.$measure" | awk '{print toupper($1)}' > $tmp/title
file=$tmp/barplot
rm -rf $file
for ii in `echo $list`
do
	i=`echo $ii | cut -f 1 -d "."`
	j=`echo $ii | cut -f 2 -d "."`
	a1=`cat ../$i.tophat/$dftA/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	a2=`cat ../$i.star/$dftA/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	b1=`cat ../$i.tophat/$dftB/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	b2=`cat ../$i.star/$dftB/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	c1=`cat ../$i.tophat/$dftC/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	c2=`cat ../$i.star/$dftC/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
	echo "$j $a1 $a2 $b1 $b2 $c1 $c2" >> $file
done
R CMD BATCH barplot.6.R
mv $tmp/pre.pdf $tmp/sim.$measure.pre.pdf
mv $tmp/sen.pdf $tmp/sim.$measure.sen.pdf
mv $file $tmp/sim.$measure.sen.summary
