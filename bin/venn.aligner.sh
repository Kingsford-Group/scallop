#!/bin/bash

file=./venn.B502.B
rm -rf $file

for i in `cat list10`
do
	echo $i
	./run.venn.aligner.sh $i >> $file
done
