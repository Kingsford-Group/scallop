#!/bin/bash

IFS="
"
for i in `cat genes.list`
do
	a=`echo "$i" | cut -f 1 -d " "`
	b=`echo "$i" | cut -f 2 -d " "`

	../scallop -i gtf/$a.gtf
done
