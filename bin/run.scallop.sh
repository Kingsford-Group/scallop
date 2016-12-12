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

{ time ./scallop -i $bam -o $dir/scallop.gtf $4 > $dir/scallop.log; } 2> $dir/time.log

./gtfcompare $gtf $dir/scallop.gtf 2 > $dir/cmp.roc

mv $dir/scallop.gtf $dir/scallop0.gtf
./gtfformat $dir/scallop0.gtf $dir/scallop.gtf

cd $dir
cuffcompare -r $gtf $dir/scallop.gtf
cd -

refallsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
refmulsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g' `
./cuffroc $dir/cuffcmp.scallop.gtf.tmap $dir/scallop.gtf $refallsize > $dir/cuffcmp.all.roc
#./cuffroc $dir/cuffcmp.scallop.gtf.tmap $dir/scallop.gtf $refmulsize 2 > $dir/cuffcmp.mul.roc
