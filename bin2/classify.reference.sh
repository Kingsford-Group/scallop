#!/bin/bash

IFS="
"
ref=./reference9
gtf=./gtf

#mkdir -p $ref
#rm -rf $ref/*.gtf

for i in `cat genes.list`
do
	a=`echo "$i" | cut -f 1 -d " "`
	b=`echo "$i" | cut -f 2 -d " "`
	c=`echo "$i" | cut -f 3 -d " "`

	if [ "$b" -ge 10 ]
	then
		continue
	fi

	echo $a $b $c

	if [ "$c" == "TRIVIAL" ]
	then
		cat $gtf/$a.gtf >> $ref/$b.trivial.gtf
	elif [ "$c" == "NORMAL" ]
	then
		cat $gtf/$a.gtf >> $ref/$b.normal.gtf
	fi
done
