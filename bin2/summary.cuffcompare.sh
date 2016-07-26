#!/bin/bash

IFS="
"

for i in `ls $1 | grep -v gtf | grep -v log | grep -v track | grep -v loci`
do
	x1=`cat $1/$i | grep mRNA | grep Query | cut -c 26-27`
	x2=`cat $1/$i | grep mRNA | grep Reference | cut -c 26-27`
	x3=`cat $1/$i | grep chain | grep Intron | cut -c 22-26`
	x4=`cat $1/$i | grep chain | grep Intron | cut -c 28-32`

	echo $i $x1 $x2 $x3 $x4
done
