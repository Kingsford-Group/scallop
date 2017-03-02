#!/bin/bash

tag="B432"
measure="all"

verA="scallop".$tag
verB="stringtie.L200.C0.0"
verC="transcomb"

dftA="scallop".$tag
dftB="stringtie.1.3.1c"
dftC="transcomb"

list="GSM981256 \
	  GSM981244 \
	  GSM984609 \
	  SRR307911 \
	  SRR387662"

tmp=./plot
mkdir -p $tmp

for i in `echo $list`
do
	echo "$i.$measure" | awk '{print toupper($i)}' > $tmp/title

	cat ../$i.tophat/$verA/cuffcmp.$measure.roc | grep ROC: > $tmp/roc1
	cat ../$i.tophat/$verB/cuffcmp.$measure.roc | grep ROC: > $tmp/roc2
	cat ../$i.tophat/$verC/cuffcmp.$measure.roc | grep ROC: > $tmp/roc3
	cat ../$i.star/$verA/cuffcmp.$measure.roc | grep ROC: > $tmp/roc4
	cat ../$i.star/$verB/cuffcmp.$measure.roc | grep ROC: > $tmp/roc5
	cat ../$i.star/$verC/cuffcmp.$measure.roc | grep ROC: > $tmp/roc6
	cat ../$i.hisat/$verA/cuffcmp.$measure.roc | grep ROC: > $tmp/roc7
	cat ../$i.hisat/$verB/cuffcmp.$measure.roc | grep ROC: > $tmp/roc8

	cat ../$i.tophat/$dftA/cuffcmp.$measure.roc | head -n 1 > $tmp/rocA
	cat ../$i.tophat/$dftB/cuffcmp.$measure.roc | head -n 1 > $tmp/rocB
	cat ../$i.tophat/$dftC/cuffcmp.$measure.roc | head -n 1 > $tmp/rocC
	cat ../$i.star/$dftA/cuffcmp.$measure.roc | head -n 1 > $tmp/rocD
	cat ../$i.star/$dftB/cuffcmp.$measure.roc | head -n 1 > $tmp/rocE
	cat ../$i.star/$dftC/cuffcmp.$measure.roc | head -n 1 > $tmp/rocF
	cat ../$i.hisat/$dftA/cuffcmp.$measure.roc | head -n 1 > $tmp/rocG
	cat ../$i.hisat/$dftB/cuffcmp.$measure.roc | head -n 1 > $tmp/rocH

	R CMD BATCH ./ROC.aligner.3.R
	mv $tmp/ROC.pdf $tmp/$i.$measure.pdf
done

list="EP \
	  ER"

for i in `echo $list`
do
	echo "$i.$measure" | awk '{print toupper($i)}' > $tmp/title

	cat ../$i.tophat/$verA/cuffcmp.$measure.roc | grep ROC: > $tmp/roc1
	cat ../$i.tophat/$verB/cuffcmp.$measure.roc | grep ROC: > $tmp/roc2
	cat ../$i.tophat/$verC/cuffcmp.$measure.roc | grep ROC: > $tmp/roc3
	cat ../$i.star/$verA/cuffcmp.$measure.roc | grep ROC: > $tmp/roc4
	cat ../$i.star/$verB/cuffcmp.$measure.roc | grep ROC: > $tmp/roc5
	cat ../$i.star/$verC/cuffcmp.$measure.roc | grep ROC: > $tmp/roc6

	cat ../$i.tophat/$dftA/cuffcmp.$measure.roc | head -n 1 > $tmp/rocA
	cat ../$i.tophat/$dftB/cuffcmp.$measure.roc | head -n 1 > $tmp/rocB
	cat ../$i.tophat/$dftC/cuffcmp.$measure.roc | head -n 1 > $tmp/rocC
	cat ../$i.star/$dftA/cuffcmp.$measure.roc | head -n 1 > $tmp/rocD
	cat ../$i.star/$dftB/cuffcmp.$measure.roc | head -n 1 > $tmp/rocE
	cat ../$i.star/$dftC/cuffcmp.$measure.roc | head -n 1 > $tmp/rocF

	R CMD BATCH ./ROC.aligner.2.R
	mv $tmp/ROC.pdf $tmp/$i.$measure.pdf
done
