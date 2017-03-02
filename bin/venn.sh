#!/bin/bash

file=./venn.W505.1.3.2d
rm -rf $file

for i in `cat list10`
do
	nohup ./run.venn.sh $i & 
done
