#!/bin/bash

if [ "$#" != 2 ]
then
	echo "usage $0 <gene ID (without scallop) of scallop> <output-dir>"
	exit
fi

file=./scallop.gtf
tmp=/tmp/scallop.$1

mkdir -p $2

cat $file | grep "gene_id \"scallop.$1\"" | grep -v abundance > $tmp

for i in `cat $tmp | cut -f 9 | cut -f 4 -d " " | cut -f 2 -d "\"" | sort | uniq`
do
	cat $tmp | grep "transcript_id \"$i\"" > $2/$i
done
