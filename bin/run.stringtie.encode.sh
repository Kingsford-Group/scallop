#!/bin/bash

if [ "$#" -ne 3 ] && [ "$#" -ne 4 ]
then
	echo "usage $0 id genome version [parameters]"
	exit
fi


encode=/home/mingfus/data/repositories/scallop/bin/ENCODE
dir=$encode/stringtie.$3/$1
bam=$encode/bam/$1.bam

mkdir -p $dir

gtf=/home/mingfus/data/repositories/scallop/bin

if [ "$2" == "GRCh38" ]; then
	gtf="$gtf"/GRCh38.gtf
elif [ "$2" == "hg19" ]; then
	gtf="$gtf"/GRCh37.gtf
elif [ "$2" == "mm10" ]; then
	gtf="$gtf"/GRCm38.gtf
elif [ "$2" == "mm9" ]; then
	gtf="$gtf"/GRCm37.gtf
fi
	
{ /usr/bin/time -v stringtie $bam -o $dir/stringtie.gtf $4 > $dir/stringtie.log; } 2> $dir/time.log

cd $dir
gffcompare -o gffmul -r $gtf $dir/stringtie.gtf -M -N
cd -
