#!/bin/bash

tag="B460"
measure="all"

dftA="scallop".$tag
dftB="stringtie.1.3.1c"
dftC="transcomb"

tmp=./plot
mkdir -p $tmp

for s in `echo "tophat star"`
do
	echo "SIMULATION.$s" | awk '{print toupper($1)}' > $tmp/title
	file=$tmp/barplot
	rm -rf $file
	for i in `echo "EP ER"`
	do
		a=`cat ../$i.$s/$dftA/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
		b=`cat ../$i.$s/$dftB/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
		c=`cat ../$i.$s/$dftC/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
		echo "$i $a $b $c" >> $file
	done
	R CMD BATCH barplot.3.R
	mv $tmp/pre.pdf $tmp/sim.$measure.$s.pre.pdf
	mv $tmp/sen.pdf $tmp/sim.$measure.$s.sen.pdf
done


list="GSM981256.ST1 \
	  GSM981244.ST2 \
	  GSM984609.ST3 \
	  SRR387662.TC1 \
	  SRR307911.TC2"

tmp=./plot
mkdir -p $tmp

for s in `echo "tophat star"`
do
	echo "REAL.DATA.$s" | awk '{print toupper($1)}' > $tmp/title
	file=$tmp/barplot
	rm -rf $file
	for ii in `echo $list`
	do
		i=`echo $ii | cut -f 1 -d "."`
		j=`echo $ii | cut -f 2 -d "."`
		a=`cat ../$i.$s/$dftA/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
		b=`cat ../$i.$s/$dftB/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
		c=`cat ../$i.$s/$dftC/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
		echo "$j $a $b $c" >> $file
	done
	R CMD BATCH barplot.3.R
	mv $tmp/pre.pdf $tmp/sra.$measure.$s.pre.pdf
	mv $tmp/sen.pdf $tmp/sra.$measure.$s.sen.pdf
done

for s in `echo "hisat"`
do
	echo "REAL.DATA.$s" | awk '{print toupper($1)}' > $tmp/title
	file=$tmp/barplot
	rm -rf $file
	for ii in `echo $list`
	do
		i=`echo $ii | cut -f 1 -d "."`
		j=`echo $ii | cut -f 2 -d "."`
		a=`cat ../$i.$s/$dftA/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
		b=`cat ../$i.$s/$dftB/cuffcmp."$measure".roc | head -n 1 | awk '{print $10,$13,$16}'`
		echo "$j $a $b" >> $file
	done
	R CMD BATCH barplot.2.R
	mv $tmp/pre.pdf $tmp/sra.$measure.$s.pre.pdf
	mv $tmp/sen.pdf $tmp/sra.$measure.$s.sen.pdf
done
