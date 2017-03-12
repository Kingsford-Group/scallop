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

#{ /usr/bin/time -v ./scallop -i $bam -o $dir/scallop.gtf $4 > $dir/scallop.log; } 2> $dir/time.log
#
#mv $dir/scallop.gtf $dir/scallop0.gtf
#./gtfformat format $dir/scallop0.gtf $dir/scallop.gtf
#
#if [ "$2" == "hisat" ]; then
#	cat $dir/scallop.gtf | sed 's/^/chr/g' > $dir/scallop.tmp.xxx.gtf
#	mv $dir/scallop.tmp.xxx.gtf $dir/scallop.gtf
#fi
#
#cd $dir
#gffcompare -o gffmul -r $gtf $dir/scallop.gtf -M -N
#cd -
#
#refmulsize=`cat $dir/gffmul.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g' `
#./gtfcuff roc $dir/gffmul.scallop.gtf.tmap $refmulsize > $dir/gffmul.roc

#./gtfcuff classify $dir/gffmul.scallop.gtf.tmap $dir/scallop.gtf > $dir/gffmul.class

id=`echo $1 | cut -f 1 -d "."`
quantfile=/home/mingfus/data/transcriptomics/SRA/"$id".all/salmon/salmon.quant/quant.sf
./gtfcuff acc-quant $dir/gffmul.scallop.gtf.tmap $quantfile 0.1 > $dir/gffmul.quant
