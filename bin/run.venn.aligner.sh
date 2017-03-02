#!/bin/bash

if [ "$#" -ne 1 ]
then
	echo "usage $0 dataset"
	exit
fi

sc="scallop.B505.0.01"
st="stringtie.1.3.2d.0.01"
tc="transcomb.0.01"

dir=`pwd`/$1.V505.1.3.2d
gtf=`pwd`/p2_sorted.gtf

roc=`pwd`/gtfcuff
commgtf=`pwd`/common.gtf.sh

mkdir -p $dir


ln -sf `pwd`/$1.tophat/$sc/scallop.gtf $dir/scallop.1.gtf
ln -sf `pwd`/$1.star/$sc/scallop.gtf $dir/scallop.2.gtf
ln -sf `pwd`/$1.hisat/$sc/scallop.gtf $dir/scallop.3.gtf

ln -sf `pwd`/$1.tophat/$st/st.gtf $dir/stringtie.1.gtf
ln -sf `pwd`/$1.star/$st/st.gtf $dir/stringtie.2.gtf
ln -sf `pwd`/$1.hisat/$st/st.gtf $dir/stringtie.3.gtf

ln -sf `pwd`/$1.tophat/$tc/TransComb.gtf $dir/transcomb.1.gtf
ln -sf `pwd`/$1.star/$tc/TransComb.gtf $dir/transcomb.2.gtf

cd $dir

for k in `echo "1 2 3"`
do
	gffcompare -r $gtf scallop.$k.gtf -M -N
	refsize=`cat $dir/gffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
	$roc split gffcmp.scallop.$k.gtf.tmap scallop.$k.gtf scallop.$k.true.gtf false.gtf
	$roc roc gffcmp.scallop.$k.gtf.tmap $refsize > scallop.$k.roc
	rm gffcmp.* false.gtf
done

for k in `echo "1 2 3"`
do
	gffcompare -r $gtf stringtie.$k.gtf -M -N
	refsize=`cat $dir/gffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
	$roc split gffcmp.stringtie.$k.gtf.tmap stringtie.$k.gtf stringtie.$k.true.gtf false.gtf
	$roc roc gffcmp.stringtie.$k.gtf.tmap $refsize > stringtie.$k.roc
	rm gffcmp.* false.gtf
done

for k in `echo "1 2"`
do
	gffcompare -r $gtf transcomb.$k.gtf -M -N
	refsize=`cat $dir/gffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
	$roc split gffcmp.transcomb.$k.gtf.tmap transcomb.$k.gtf transcomb.$k.true.gtf false.gtf
	$roc roc gffcmp.transcomb.$k.gtf.tmap $refsize > transcomb.$k.roc
	rm gffcmp.* false.gtf
done

$commgtf $dir/scallop.1.true.gtf $dir/scallop.2.true.gtf > $dir/scallop.12.true.gtf
$commgtf $dir/scallop.1.true.gtf $dir/scallop.3.true.gtf > $dir/scallop.13.true.gtf
$commgtf $dir/scallop.2.true.gtf $dir/scallop.3.true.gtf > $dir/scallop.23.true.gtf
$commgtf $dir/scallop.1.true.gtf $dir/scallop.23.true.gtf > $dir/scallop.123.true.gtf

$commgtf $dir/stringtie.1.true.gtf $dir/stringtie.2.true.gtf > $dir/stringtie.12.true.gtf
$commgtf $dir/stringtie.1.true.gtf $dir/stringtie.3.true.gtf > $dir/stringtie.13.true.gtf
$commgtf $dir/stringtie.2.true.gtf $dir/stringtie.3.true.gtf > $dir/stringtie.23.true.gtf
$commgtf $dir/stringtie.1.true.gtf $dir/stringtie.23.true.gtf > $dir/stringtie.123.true.gtf

$commgtf $dir/transcomb.1.true.gtf $dir/transcomb.2.true.gtf > $dir/transcomb.12.true.gtf

xx=$1
for k in `echo "1 2 3 12 13 23 123"`
do
	x=`cat scallop.$k.true.gtf | awk '$3 == "transcript"' | wc -l`
	xx="$xx $x"
done

for k in `echo "1 2 3 12 13 23 123"`
do
	x=`cat stringtie.$k.true.gtf | awk '$3 == "transcript"' | wc -l`
	xx="$xx $x"
done

for k in `echo "1 2 12"`
do
	x=`cat transcomb.$k.true.gtf | awk '$3 == "transcript"' | wc -l`
	xx="$xx $x"
done
echo $xx

cd -
