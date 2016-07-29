#!/bin/bash

ref=./reference9
sc=./scallop9a/
st=./stringtie9a/

for i in `seq 1 9`
do
	./collect.gtf.sh $sc/$i.normal full > $sc/$i.normal.gtf
	./collect.gtf.sh $sc/$i.trivial full > $sc/$i.trivial.gtf
	a1=`../gtfcompare $ref/$i.normal.gtf $sc/$i.normal.gtf 2 | cut -f 3,6,9,12,15 -d " "`
	a2=`../gtfcompare $ref/$i.trivial.gtf $sc/$i.trivial.gtf 2 | cut -f 3,6,9,12,15 -d " "`

	./collect.gtf.sh $st/$i.normal STRG > $st/$i.normal.gtf
	./collect.gtf.sh $st/$i.trivial STRG > $st/$i.trivial.gtf
	b1=`../gtfcompare $ref/$i.normal.gtf $st/$i.normal.gtf 2 | cut -f 3,6,9,12,15 -d " "`
	b2=`../gtfcompare $ref/$i.trivial.gtf $st/$i.trivial.gtf 2 | cut -f 3,6,9,12,15 -d " "`

	echo "$i TRIVIAL $a1 $b1" | awk '{print $1, $2, $3, $4, $9, $5, $10, $6, $11, $7, $12}'
	echo "$i NORMAL $a2 $b2" | awk '{print $1, $2, $3, $4, $9, $5, $10, $6, $11, $7, $12}'
done
