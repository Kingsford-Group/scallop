#!/bin/bash

if [ "$#" -ne 3 ] && [ "$#" -ne 4 ]
then
	echo "usage $0 dataset aligner version [parameters]"
	exit
fi

dir=`pwd`/$1/stringtie."$3"
mkdir -p $dir

bam=`pwd`/$1/$2/"$2".sort.bam
gtf=`pwd`/$1/expression.gtf

{ time stringtie $bam -o $dir/st.gtf $4 > $dir/st.log; } 2> $dir/time.log

./gtfcompare $gtf $dir/st.gtf 2 > $dir/cmp.roc

mv $dir/st.gtf $dir/st0.gtf
./gtfformat $dir/st0.gtf $dir/st.gtf

cd $dir
cuffcompare -r $gtf $dir/st.gtf
cd -

refallsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
refmulsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g' `
./cuffroc $dir/cuffcmp.st.gtf.tmap $dir/st.gtf $refallsize 1 > $dir/cuffcmp.all.roc
./cuffroc $dir/cuffcmp.st.gtf.tmap $dir/st.gtf $refmulsize 2 > $dir/cuffcmp.mul.roc
