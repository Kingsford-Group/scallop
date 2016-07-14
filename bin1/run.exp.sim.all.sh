#!/bin/bash

dir=`pwd`/
list="mouse zebrafish"

for i in `echo $list`
do
	echo $i
	./run.exp.sim.sh $i
done
