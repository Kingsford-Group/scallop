#!/bin/bash

if [ "$#" != 1 ]
then
	echo "usage $0 <gene-name>"
	exit
fi

mkdir -p tex
ln -sf ../llncs.cls tex/

./scallop -i ./expression.gtf -g $1 -a full -o $1.x.gtf -t > $1.x.log 

mkdir -p tex
rm -rf tex/$1.*.tex
mv $1.*.tex tex
./create.pdf.sh $1


./scallop -i ./genes/$1.sort.bam -a full -o $1.y.gtf -t > $1.y.log 
rm -rf tex/bundle.0.*.tex
mv bundle.0.*.tex tex
./create.pdf.sh bundle.0
