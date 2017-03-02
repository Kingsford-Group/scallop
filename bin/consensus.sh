#!/bin/bash

for i in `cat list10`
do
	echo $i
	./run.consensus.sh $i
done
