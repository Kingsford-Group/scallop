#!/bin/bash

tag="B432"
verA="scallop".$tag
verB="stringtie.L200.C0.0"
verC="transcomb"

dftA="scallop".$tag
dftB="stringtie.1.3.1c"
dftC="transcomb"

list="EP.tophat \
	  EP.star \
	  ER.tophat \
	  ER.star \
	  GSM981256.tophat \
	  GSM981256.star \
	  GSM981244.tophat \
	  GSM981244.star \
	  GSM984609.tophat \
	  GSM984609.star \
	  SRR387662.tophat \
	  SRR387662.star \
	  SRR307911.tophat \
	  SRR307911.star"

tmp=./plot
mkdir -p $tmp

for i in `echo $list`
do
	echo $i | awk '{print toupper($1)}' > $tmp/title
	cat ../$i/$verA/cuffcmp.all.roc | grep ROC: > $tmp/roc1
	cat ../$i/$verB/cuffcmp.all.roc | grep ROC: > $tmp/roc2
	cat ../$i/$verC/cuffcmp.all.roc | grep ROC: > $tmp/roc3
	cat ../$i/$verA/cuffcmp.mul.roc | grep ROC: > $tmp/roc4
	cat ../$i/$verB/cuffcmp.mul.roc | grep ROC: > $tmp/roc5
	cat ../$i/$verC/cuffcmp.mul.roc | grep ROC: > $tmp/roc6

	cat ../$i/$dftA/cuffcmp.all.roc | head -n 1 > $tmp/rocA
	cat ../$i/$dftB/cuffcmp.all.roc | head -n 1 > $tmp/rocB
	cat ../$i/$dftC/cuffcmp.all.roc | head -n 1 > $tmp/rocC
	cat ../$i/$dftA/cuffcmp.mul.roc | head -n 1 > $tmp/rocD
	cat ../$i/$dftB/cuffcmp.mul.roc | head -n 1 > $tmp/rocE
	cat ../$i/$dftC/cuffcmp.mul.roc | head -n 1 > $tmp/rocF

	R CMD BATCH ./ROC.all.mul.3.R
	mv $tmp/ROC.pdf $tmp/$i.pdf
done

list="GSM981256.hisat \
	  GSM981244.hisat \
	  GSM984609.hisat \
	  SRR387662.hisat \
	  SRR307911.hisat"

for i in `echo $list`
do
	echo $i | awk '{print toupper($1)}' > $tmp/title

	cat ../$i/$verA/cuffcmp.all.roc | grep ROC: > $tmp/roc1
	cat ../$i/$verB/cuffcmp.all.roc | grep ROC: > $tmp/roc2
	cat ../$i/$verA/cuffcmp.mul.roc | grep ROC: > $tmp/roc4
	cat ../$i/$verB/cuffcmp.mul.roc | grep ROC: > $tmp/roc5

	cat ../$i/$dftA/cuffcmp.all.roc | head -n 1 > $tmp/rocA
	cat ../$i/$dftB/cuffcmp.all.roc | head -n 1 > $tmp/rocB
	cat ../$i/$dftA/cuffcmp.mul.roc | head -n 1 > $tmp/rocD
	cat ../$i/$dftB/cuffcmp.mul.roc | head -n 1 > $tmp/rocE

	R CMD BATCH ./ROC.all.mul.2.R
	mv $tmp/ROC.pdf $tmp/$i.pdf
done
