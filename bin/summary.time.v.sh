#!/bin/bash

tag=B502.B
ver=`stringtie --version`

list="GSM981256.ST1 \
	  GSM981244.ST2 \
	  GSM984609.ST3 \
	  SRR387661.TC1 \
	  SRR307911.TC2 \
	  SRR545723.SC1 \
	  SRR315323.SC2 \
	  SRR307903.SC3 \
	  SRR315334.SC4 \
	  SRR534307.SC5"

for j in `echo $list`
do
	i=`echo $j | cut -f 1 -d "."`
	d=`echo $j | cut -f 2 -d "."`

	k=""
	for m in `echo "scallop.B505.A.0.01 stringtie.1.3.2d.A.0.01"`
	do
		for a in `echo "tophat star hisat"`
		do
			u=`cat $i.$a/$m/time.log | grep time | grep User | awk '{print $4}'`
			s=`cat $i.$a/$m/time.log | grep time | grep Sys  | awk '{print $4}'`
			t=`echo "$u $s" | awk '{print $1 + $2}'`
			k="$k $t"
		done
	done

	for m in `echo "transcomb.A."`
	do
		for a in `echo "tophat star"`
		do
			u=`cat $i.$a/$m/time.log | grep time | grep User | awk '{print $4}'`
			s=`cat $i.$a/$m/time.log | grep time | grep Sys  | awk '{print $4}'`
			t=`echo "$u $s" | awk '{print $1 + $2}'`
			k="$k $t"
		done
	done

	echo $d $k
done
