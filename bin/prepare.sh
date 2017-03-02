#!/bin/bash

dir=`pwd`/
sra=/home/mingfus/data/transcriptomics/SRA

for i in `cat list`
do
	echo $i

	cur=$dir/$i.tophat
	mkdir -p $cur
	cd $cur
	ln -sf $sra/$i/tophat .
	ln -sf ../human.p2_sorted.gtf expression.gtf
	cd -

	cur=$dir/$i.star
	mkdir -p $cur
	cd $cur
	ln -sf $sra/$i/star .
	ln -sf ../human.p2_sorted.gtf expression.gtf
	cd -

	cur=$dir/$i.hisat
	mkdir -p $cur
	cd $cur
	ln -sf $sra/$i/hisat .
	ln -sf ../human.p2_sorted_no.gtf expression.gtf
	cd -

done
