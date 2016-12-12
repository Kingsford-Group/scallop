#!/bin/bash

ref=./reference9
sc=./scallop9q
st=./stringtie09

for i in `seq 1 9`
do
	./collect.gtf.sh $sc/s1.$i core > $sc/s1.$i.gtf
	./collect.gtf.sh $sc/s2.$i full > $sc/s2.$i.gtf
	a1=`../gtfcompare $ref/$i.gtf $sc/s1.$i.gtf 2 | cut -f 3,6,9,12,15 -d " "`
	a2=`../gtfcompare $ref/$i.gtf $sc/s2.$i.gtf 2 | cut -f 3,6,9,12,15 -d " "`

#./collect.gtf.sh $st/$i STRG > $st/$i.gtf
	b1=`../gtfcompare $ref/$i.gtf $st/$i.gtf 2 | cut -f 3,6,9,12,15 -d " "`

	echo "$i ALL $a1 $a2 $b1" | awk '{print $1, $2, $3, $4, $9, $14, $5, $10, $15, $6, $11, $16, $7, $12, $17}'
done
