#!/bin/bash

if [ "$#" != 1 ]
then
	echo "usage $0 <gene-name>"
	exit
fi

mkdir -p tex
ln -sf ../llncs.cls tex/

ln -sf ../bam/$1.bam .
samtools index $1.bam

stringtie $1.bam -o $1.st.gtf
cuffcompare -r ../gtf/$1.gtf $1.st.gtf

../../scallop -i $1.bam -a full -o $1.s1.gtf -t > $1.s1.log
mkdir -p tex
rm -rf tex/bundle.0.*.tex
mv bundle.0.*.tex tex
./create.pdf.sh bundle.0

../../scallop -i ../gtf/$1.gtf -a full -o $1.s2.gtf -t > $1.s2.log
mkdir -p tex
rm -rf tex/$1.*.tex
mv $1.*.tex tex
./create.pdf.sh $1
