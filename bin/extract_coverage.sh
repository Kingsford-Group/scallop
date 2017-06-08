#!/bin/bash

sc=scallop23
cov5=./coverage5
cov3=./coverage3
tmp1=/tmp/tmp1
tmp2=/tmp/tmp2

for i in `cat list23`
do
	cat $sc/$i.log | grep 5end | cut -f 4 -d " " > $tmp1
	l=`cat $tmp1 | wc -l`

	echo $i $l

	if [ "$l" == "500" ]
	then
		paste $cov5 $tmp1 > $tmp2
		mv $tmp2 $cov5
	fi

	cat $sc/$i.log | grep 3end | cut -f 4 -d " " > $tmp1
	l=`cat $tmp1 | wc -l`

	echo $i $l

	if [ "$l" == "500" ]
	then
		paste $cov3 $tmp1 > $tmp2
		mv $tmp2 $cov3
	fi
done
