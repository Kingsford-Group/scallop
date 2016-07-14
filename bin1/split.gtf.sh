#!/bin/bash

if [ "$#" != 3 ]
then
	echo "usage $0 <gene ID> <output-dir> <gtf.file>"
	exit
fi

file=$3
tmp=/tmp/scallop.$1

mkdir -p $2

cat $file | grep "gene_id \"$1\"" | grep -v abundance | grep -v FPKM > $tmp

for i in `cat $tmp | cut -f 9 | cut -f 4 -d " " | cut -f 2 -d "\"" | sort | uniq`
do
	cat $tmp | grep "transcript_id \"$i\"" > $2/"$i".gtf
done
