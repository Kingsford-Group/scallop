#!/bin/bash

if [ "$#" -ne 2 ]
then
	echo "usage $0 dataset aligner"
	exit
fi

ver=`stringtie --version`
dir=`pwd`/$1/stringtie."$ver"
mkdir -p $dir


bam=`pwd`/$1/$2/"$2".sort.bam
gtf=`pwd`/$1/expression.gtf

#{ time stringtie $bam -o $dir/st.gtf > $dir/st.log; } 2> $dir/time.log

./gtfcompare $gtf $dir/st.gtf 2 > $dir/cmp1.log
./gtfcompare $dir/st.gtf $gtf 2 > $dir/cmp2.log

cd $dir
cuffcompare -r $gtf $dir/st.gtf
cd -
