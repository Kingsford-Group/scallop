#!/bin/bash

if [ "$#" -ne 1 ]
then
	echo "usage $0 dataset"
	exit
fi

id="B502.C"
dir=`pwd`/$1.C502.C
gtf=`pwd`/p2_sorted.gtf
roc=`pwd`/cuffroc

mkdir -p $dir

sc="scallop.$id"
st="stringtie.1.3.1c"
tc="transcomb"

ln -sf `pwd`/$1.tophat/$sc/scallop.gtf $dir/scallop.1.gtf
ln -sf `pwd`/$1.star/$sc/scallop.gtf $dir/scallop.2.gtf
ln -sf `pwd`/$1.hisat/$sc/scallop.gtf $dir/scallop.3.gtf

ln -sf `pwd`/$1.tophat/$st/st.gtf $dir/stringtie.1.gtf
ln -sf `pwd`/$1.star/$st/st.gtf $dir/stringtie.2.gtf
ln -sf `pwd`/$1.hisat/$st/st.gtf $dir/stringtie.3.gtf

ln -sf `pwd`/$1.tophat/$tc/TransComb.gtf $dir/transcomb.1.gtf
ln -sf `pwd`/$1.star/$tc/TransComb.gtf $dir/transcomb.2.gtf

./common.gtf.sh $dir/scallop.1.gtf $dir/scallop.2.gtf > $dir/scallop.12.gtf
./common.gtf.sh $dir/scallop.1.gtf $dir/scallop.3.gtf > $dir/scallop.13.gtf
./common.gtf.sh $dir/scallop.2.gtf $dir/scallop.3.gtf > $dir/scallop.23.gtf
./common.gtf.sh $dir/scallop.1.gtf $dir/scallop.23.gtf > $dir/scallop.123.gtf

./common.gtf.sh $dir/stringtie.1.gtf $dir/stringtie.2.gtf > $dir/stringtie.12.gtf
./common.gtf.sh $dir/stringtie.1.gtf $dir/stringtie.3.gtf > $dir/stringtie.13.gtf
./common.gtf.sh $dir/stringtie.2.gtf $dir/stringtie.3.gtf > $dir/stringtie.23.gtf
./common.gtf.sh $dir/stringtie.1.gtf $dir/stringtie.23.gtf > $dir/stringtie.123.gtf

./common.gtf.sh $dir/transcomb.1.gtf $dir/transcomb.2.gtf > $dir/transcomb.12.gtf

cd $dir

for k in `echo "12 13 23 123"`
do
	cuffcompare -r $gtf scallop.$k.gtf
	refsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
	$roc cuffcmp.scallop.$k.gtf.tmap scallop.$k.gtf $refsize > scallop.$k.roc
	rm cuffcmp*
done

for k in `echo "12 13 23 123"`
do
	cuffcompare -r $gtf stringtie.$k.gtf
	refsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
	$roc cuffcmp.stringtie.$k.gtf.tmap stringtie.$k.gtf $refsize > stringtie.$k.roc
	rm cuffcmp*
done

for k in `echo "12"`
do
	cuffcompare -r $gtf transcomb.$k.gtf
	refsize=`cat $dir/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
	$roc cuffcmp.transcomb.$k.gtf.tmap transcomb.$k.gtf $refsize > transcomb.$k.roc
	rm cuffcmp*
done

cd -
