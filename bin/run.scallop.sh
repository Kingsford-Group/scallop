#!/bin/bash

if [ "$#" -ne 3 ] && [ "$#" -ne 4 ]
then
	echo "usage $0 dataset aligner version [parameters]"
	exit
fi

dir=`pwd`/$1/scallop.$3
mkdir -p $dir

bam=`pwd`/$1/$2/"$2".sort.bam
gtf=`pwd`/$1/expression.gtf

{ /usr/bin/time -v ./scallop -i $bam -o $dir/scallop.gtf $4 > $dir/scallop.log; } 2> $dir/time.log

mv $dir/scallop.gtf $dir/scallop0.gtf
./gtfformat $dir/scallop0.gtf $dir/scallop.gtf

cd $dir
gffcompare -o gffmul -r $gtf $dir/scallop.gtf -M -N
cd -

refmulsize=`cat $dir/gffcmp.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g' `
./gtfcuff roc $dir/cuffcmp.scallop.gtf.tmap $refmulsize > $dir/gffmul.roc
