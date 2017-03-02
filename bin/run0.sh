#!/bin/bash

tag=$1

list="GSM981256.tophat"

for i in `echo $list`
do
	aligner=`echo $i | cut -f 2 -d "."`
	echo $i $aligner
	nohup ./run.scallop.sh $i $aligner B412 &
	nohup ./run.stringtie.sh $i $aligner &
	nohup ./run.transcomb.sh $i $aligner &
done
