#!/bin/bash

dir=`pwd`/
list="fruitfly chimpanzee mouse rat zebrafish"

for i in `echo $list`
do
	echo $i
	./run.exp.sim.2.sh $i
done
