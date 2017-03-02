#!/bin/bash

if [ "$#" -ne 2 ]
then
	echo "usage $0 <gtf1> <gtf2>"
	exit
fi

roc=/home/mingfus/data/repositories/scallop/bin/gtfcuff
dir=`pwd`/

tmp=`mktemp -d -p .`

echo $tmp

cd $dir/$tmp

ln -sf $1 aaa.gtf
ln -sf $2 bbb.gtf

gffcompare -r aaa.gtf bbb.gtf 1> /dev/null 2> /dev/null
$roc split gffcmp.bbb.gtf.tmap bbb.gtf true.gtf false.gtf 1> /dev/null 2> /dev/null
rm -rf gffcmp* false.gtf
cat true.gtf

rm -rf $dir/$tmp
