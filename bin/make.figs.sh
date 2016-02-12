#!/bin/bash

figs=./figs
mkdir -p $figs

./scallop config sort.bam > log

rm -rf ./$figs/*
mv ./sgraph*.tex $figs

cp llncs.cls $figs

cd $figs

for i in `seq 0 9`
do
	echo $i
	myepstool.sh sgraph"$i"a
	myepstool.sh sgraph"$i"b
	rm -rf sgraph"$i"*.aux sgraph"$i"*.eps sgraph"$i"*.dvi sgraph"$i"*.log
done
