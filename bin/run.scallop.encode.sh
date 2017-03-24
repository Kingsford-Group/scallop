#!/bin/bash

if [ "$#" -ne 3 ] && [ "$#" -ne 4 ]
then
	echo "usage $0 id genome version [parameters]"
	exit
fi


encode=/home/mingfus/data/repositories/scallop/bin/ENCODE
dir=$encode/scallop.$3/$1
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
	
{ /usr/bin/time -v ./scallop -i $bam -o $dir/scallop.gtf $4 > $dir/scallop.log; } 2> $dir/time.log
mv $dir/scallop.gtf $dir/scallop0.gtf
./gtfformat format $dir/scallop0.gtf $dir/scallop.gtf

cd $dir
gffcompare -o gffmul -r $gtf $dir/scallop.gtf -M -N
gtfcompare $gtf $dir/scallop.gtf > $dir/gtfcmp
cd -

#refmulsize=`cat $dir/gffmul.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g' `
#./gtfcuff roc $dir/gffmul.scallop.gtf.tmap $refmulsize > $dir/gffmul.roc
