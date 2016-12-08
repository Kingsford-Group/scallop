#!/bin/bash

if [ "$#" -ne 3 ] && [ "$#" -ne 4 ] && [ "$#" -ne 2 ]
then
	echo "usage $0 dataset aligner [version] [parameters]"
	exit
fi

dir=`pwd`/$1/transcomb

if [ "$#" -ge 3 ]
then
	dir=`pwd`/$1/transcomb.$3
fi

mkdir -p $dir

bam=`pwd`/$1/$2/"$2".sort.bam
gtf=`pwd`/$1/expression.gtf

#cd $dir
#{ time TransComb -b $bam -o $dir $4 > $dir/transcomb.log; } 2> $dir/time.log
#cd -
#
#./gtfcompare $gtf $dir/TransComb.gtf 2 > $dir/cmp.roc
#
#mv $dir/TransComb.gtf $dir/TransComb0.gtf
#./gtfformat $dir/TransComb0.gtf $dir/TransComb.gtf
#
#cd $dir
#cuffcompare -r $gtf $dir/TransComb.gtf
#cd -

refallsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
refmulsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g'`
./cuffroc $dir/cuffcmp.TransComb.gtf.tmap $dir/TransComb.gtf $refallsize 1 > $dir/cuffcmp.all.roc
./cuffroc $dir/cuffcmp.TransComb.gtf.tmap $dir/TransComb.gtf $refmulsize 2 > $dir/cuffcmp.mul.roc
