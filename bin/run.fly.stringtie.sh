#!/bin/bash

if [ "$#" -ne 1 ]; then
	echo "usage $0 ID"
	exit
fi

fly=/home/mingfus/data/transcriptomics/fly/
gtf=$fly/BDGP6.85.gtf
dir=`pwd`/fly/$1/stringtie
mkdir -p $dir

bam=$fly/hisat/$1/hisat.sort.bam

{ /usr/bin/time -v stringtie0 $bam -o $dir/stringtie.gtf > $dir/stringtie.log; } 2> $dir/time.log

mv $dir/stringtie.gtf $dir/stringtie0.gtf
./gtfformat format $dir/stringtie0.gtf $dir/stringtie.gtf

cd $dir
gffcompare -o gffmul -r $gtf $dir/stringtie.gtf -M -N
cd -

refmulsize=`cat $dir/gffmul.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g' `
./gtfcuff roc $dir/gffmul.stringtie.gtf.tmap $refmulsize > $dir/gffmul.roc
