#!/bin/bash

if [ "$#" -ne 1 ]
then
	echo "usage: $0 cuffcompare.stats"
	exit
fi

total1=`cat $1 | grep mRNAs | grep Query | awk '{print $9}' | sed 's/(//g'`
total2=`cat $1 | grep mRNAs | grep Refer | awk '{print $9}' | sed 's/(//g'`
correct=`cat $1 | grep Matching | grep intron | awk '{print $4}'`

if [ "$total1" == "" ]
then
	total1="1"
fi
if [ "$total2" == "" ]
then
	total2="1"
fi
if [ "$correct" == "" ]
then
	correct="0"
fi

echo "$total1 $total2 $correct" | awk '{printf("%d %d %.2lf %.2lf\n", $1, $3, 100*$3/$2, 100*$3/$1)}'
