#!/bin/bash

for i in `cat list10`
do
	echo $i
	nohup ./run.venn.single.sh $i &
done
