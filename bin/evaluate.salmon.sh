#!/bin/bash

if [ "$#" -ne 3 ]
then
	echo "usage $0 dataset aligner version"
	exit
fi

gtf=$1/$1/salmon.gtf

dir=`pwd`/$1/scallop.$3

cd $dir
cuffcompare -r $gtf $dir/scallop.gtf
cd -

refallsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
./cuffroc $dir/cuffcmp.scallop.gtf.tmap $dir/scallop.gtf $refallsize > $dir/salmon.all.roc
