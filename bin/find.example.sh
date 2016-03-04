#!/bin/bash

dir=./examples
log=./find.log

mkdir -p $dir

for i in `seq 1 999999`
do
	echo $i

	./scallop 5 10 > $log
	myEpsTool.sh sgraph 2>/dev/null 1>/dev/null

	a=`cat $log | grep EXAMPLE`
	if [ "$a" = "EXAMPLE" ]
	then
		mv $log $dir/"$i".log
		mv sgraph.gr $dir/"$i".gr
	fi

	sleep 1
done
