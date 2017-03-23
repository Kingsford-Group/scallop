#!/bin/bash

if [ "$#" -ne 3 ] && [ "$#" -ne 4 ]
then
	echo "usage $0 dataset aligner version [parameters]"
	exit
fi

dir=`pwd`/$1/stringtie.$3
mkdir -p $dir

gtf=/home/mingfus/data/repositories/scallop/bin/gencode.v25.gtf

cd $dir
gffcompare -o gffmul -r $gtf $dir/st.gtf -M -N
cd -

refmulsize=`cat $dir/gffmul.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g' `
./gtfcuff roc $dir/gffmul.st.gtf.tmap $refmulsize > $dir/gffmul.roc

#id=`echo $1 | cut -f 1 -d "."`
#quantfile=/home/mingfus/data/transcriptomics/SRA/"$id".all/salmon/salmon.quant/quant.sf
#./gtfcuff acc-quant $dir/gffmul.st.gtf.tmap $quantfile 0.1 > $dir/gffmul.quant
