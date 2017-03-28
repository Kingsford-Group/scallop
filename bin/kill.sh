#!/bin/bash

for i in `seq 1 100`
do
	killall stringtie scallop gtfcompare gffcompare
	sleep 1
done
