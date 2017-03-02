#!/bin/bash

if [ "$#" -ne 2 ]
then
	echo "usage $0 dataset aligner"
	exit
fi

dir=`pwd`/$1/transcomb1
mkdir -p $dir

bam=`pwd`/$1/$2/"$2".sort.bam
gtf=`pwd`/$1/expression.gtf

cd $dir
#{ time TransComb -b $bam -o $dir -s first > $dir/transcomb.log; } 2> $dir/time.log
cd -

./gtfcompare $gtf $dir/TransComb.gtf 2 > $dir/cmp1.log
./gtfcompare $dir/TransComb.gtf $gtf 2 > $dir/cmp2.log

cd $dir
cuffcompare -r $gtf $dir/TransComb.gtf
cd -

exit

dir=`pwd`/$1/transcomb2
mkdir -p $dir

bam=`pwd`/$1/$2/"$2".sort.bam
gtf=`pwd`/$1/expression.gtf

cd $dir
{ time TransComb -b $bam -o $dir -s second > $dir/transcomb.log; } 2> $dir/time.log
cd -

./gtfcompare $gtf $dir/TransComb.gtf 2 > $dir/cmp1.log
./gtfcompare $dir/TransComb.gtf $gtf 2 > $dir/cmp2.log

cd $dir
cuffcompare -r $gtf $dir/TransComb.gtf
cd -
