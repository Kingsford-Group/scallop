#!/bin/bash

if [ "$#" != 2 ]
then
	echo "usage $0 <gene ID (without STRG) of stringtie> <output-dir>"
	exit
fi

file=./stringtie.gtf
tmp=/tmp/stringtie.$1

mkdir -p $2

cat $file | grep "gene_id \"STRG.$1\"" | grep -v FPKM > $tmp

for i in `cat $tmp | cut -f 9 | cut -f 4 -d " " | cut -f 2 -d "\"" | sort | uniq`
do
	cat $tmp | grep "transcript_id \"$i\"" > $2/$i
done
