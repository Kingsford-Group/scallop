#!/bin/bash

tag="P3"

list="GSM981256.tophat \
	  GSM981244.tophat \
	  GSM984609.tophat \
	  SRR307911.tophat \
	  SRR387661.tophat"

for i in `echo $list`
do
	nohup ./run.perfect.sh $i $tag &
done
