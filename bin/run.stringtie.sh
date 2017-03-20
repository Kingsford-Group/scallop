#!/bin/bash

if [ "$#" -ne 3 ] && [ "$#" -ne 4 ]
then
	echo "usage $0 dataset aligner version [parameters]"
	exit
fi

dir=`pwd`/$1/stringtie.$3
mkdir -p $dir

id=`echo $1 | cut -f 1 -d "."`
quantfile=/home/mingfus/data/transcriptomics/SRA/"$id".all/salmon/salmon.quant/quant.sf
./gtfcuff acc-quant $dir/gffmul.st.gtf.tmap $quantfile 0.1 > $dir/gffmul.quant
