#!/bin/bash

if [ "$#" -ne 3 ] && [ "$#" -ne 4 ]
then
	echo "usage $0 dataset aligner version [parameters]"
	exit
fi

dir=`pwd`/$1/transcomb.$3
mkdir -p $dir

gtf=/home/mingfus/data/repositories/scallop/bin/GRCh38.gtf
#gtf=/home/mingfus/data/repositories/scallop/bin/gencode.v25.gtf

#cd $dir
#gffcompare -o gffmul -r $gtf $dir/TransComb.gtf -M -N
#cd -

id=`echo $1 | cut -f 1 -d "."`
quantfile=/home/mingfus/data/transcriptomics/SRA/"$id".all/salmon/salmon.quant/quant.sf
./gtfcuff acc-quant $dir/gffmul.TransComb.gtf.tmap $quantfile 0.1 > $dir/gffmul.quant

#./gtfcuff classify $dir/gffmul.TransComb.gtf.tmap $dir/TransComb.gtf > $dir/gffmul.class
