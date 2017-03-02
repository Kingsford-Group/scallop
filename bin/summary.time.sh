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

pars="0.01 1 2.5 5 10 25 50 75 100"

for j in `echo $list`
do
	i=`echo $j | cut -f 1 -d "."`
	d=`echo $j | cut -f 2 -d "."`

	for m in `echo "scallop.B505 stringtie.1.3.2d transcomb"`
	do
		for a in `echo "tophat star hisat"`
		do
			if [ "$m" == "transcomb" ] && [ "$a" == "hisat" ]
			then
				continue;
			fi

			k=""
			for p in `echo $pars`
			do
				u1=`cat $i.$a/$m.$p/time.log | grep user | awk '{print $2}' | cut -f 1 -d "m"`
				u2=`cat $i.$a/$m.$p/time.log | grep user | awk '{print $2}' | cut -f 2 -d "m" | cut -f 1 -d "s"`
				s1=`cat $i.$a/$m.$p/time.log | grep sys  | awk '{print $2}' | cut -f 1 -d "m"`
				s2=`cat $i.$a/$m.$p/time.log | grep sys  | awk '{print $2}' | cut -f 2 -d "m" | cut -f 1 -d "s"`
				t=`echo "$u1 $u2 $s1 $s2" | awk '{print $1 * 60 + $2 + $3 * 60 + $4}'`
				k="$k $t"
			done
			echo "$d $m $a $k"
		done
	done
done
