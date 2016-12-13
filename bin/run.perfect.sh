#!/bin/bash

if [ "$#" -ne 2 ]
then
	echo "usage $0 dataset id"
	exit
fi

dir1=`pwd`/$1/scallop.$2
dir2=`pwd`/$1/stringtie.1.3.1c
dir3=`pwd`/$1/transcomb
gtf=`pwd`/$1/expression.gtf

dir=`pwd`/$1/perfect.$2
mkdir -p $dir

cd $dir
ln -sf $dir1/scallop.gtf .
ln -sf $dir2/st.gtf .
ln -sf $dir2/TransComb.gtf tc.gtf

cuffcompare -r $gtf scallop.gtf
refsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
../../cuffroc cuffcmp.scallop.gtf.tmap scallop.gtf $refsize > cuffcmp0.roc 


cuffcompare -r st.gtf scallop.gtf
refsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
../../cuffroc cuffcmp.scallop.gtf.tmap scallop.gtf $refsize > cuffcmp0.st.roc
mv true.gtf scallop1.gtf

cuffcompare -r $gtf scallop1.gtf 
refsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
../../cuffroc cuffcmp.scallop1.gtf.tmap scallop1.gtf $refsize > cuffcmp1.roc


cuffcompare -r tc.gtf scallop1.gtf
refsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
../../cuffroc cuffcmp.scallop1.gtf.tmap scallop1.gtf $refsize > cuffcmp1.tc.roc
mv true.gtf scallop2.gtf

cuffcompare -r $gtf scallop2.gtf 
refsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
../../cuffroc cuffcmp.scallop2.gtf.tmap scallop2.gtf $refsize > cuffcmp2.roc

cd -

