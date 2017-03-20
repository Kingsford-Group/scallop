#!/bin/bash

sc=scallop09
st=stringtie09
cp=compare09

tmp1=/tmp/tmp1
tmp2=/tmp/tmp2

for i in `seq 1 9`
do
	echo $i
	cat $sc/summary.$i | sort -k1,1 > $tmp1
	cat $st/summary.$i | sort -k1,1 > $tmp2
	paste $tmp1 $tmp2 | cut -f 1,2,3,4,5,7,8,9,10 | sed 's/-nan/0/g' | sed 's/-/0/g' | column -t > $cp/$i
done
