#!/bin/bash

IFS="
"
ref=./reference09
gtf=./gtf

mkdir -p $ref
rm -rf $ref/*.gtf

for i in `cat genes.list`
do
	a=`echo "$i" | cut -f 1 -d " "`
	b=`echo "$i" | cut -f 2 -d " "`

	if [ "$b" -ge 10 ]
	then
		continue
	fi

	echo $a $b

	cat $gtf/$a.gtf >> $ref/$b.gtf
done
