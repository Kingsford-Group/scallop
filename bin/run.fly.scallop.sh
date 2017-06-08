#!/bin/bash

if [ "$#" -ne 1 ]; then
	echo "usage $0 ID"
	exit
fi

fly=/home/mingfus/data/transcriptomics/fly/
gtf=$fly/BDGP6.85.gtf
dir=`pwd`/fly/$1/scallop
mkdir -p $dir

bam=$fly/hisat/$1/hisat.sort.bam

{ /usr/bin/time -v ./scallop -i $bam -o $dir/scallop.gtf > $dir/scallop.log; } 2> $dir/time.log

mv $dir/scallop.gtf $dir/scallop0.gtf
./gtfformat format $dir/scallop0.gtf $dir/scallop.gtf

cd $dir
gffcompare -o gffmul -r $gtf $dir/scallop.gtf -M -N
cd -

refmulsize=`cat $dir/gffmul.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g' `
./gtfcuff roc $dir/gffmul.scallop.gtf.tmap $refmulsize > $dir/gffmul.roc
