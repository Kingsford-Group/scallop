#!/bin/bash

file=./venn.W505.1.3.2d
rm -rf $file

for i in `cat list10`
do
	echo $i
	nohup ./run.venn.algo.sh $i & 
done
