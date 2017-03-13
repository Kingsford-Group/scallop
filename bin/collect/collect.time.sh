#!/bin/bash

if [ "$#" != "0" ]; then
	echo "usage $0"
	exit
fi

scallop="scallop.B668"
stringtie="stringtie.1.3.2d"
transcomb="transcomb"

list=../list10

for id in `cat $list`
do
	# scallop
	for aa in `echo "tophat star hisat"`
	do
		cc=""
		for abd in `echo "0.01 1 2.5 5 7.5 10 25 50 75 100"`
		do
			file="../$id.$aa/$scallop.$abd/time.log"
			t1=`cat $file | grep "User time" | awk '{print $4}'`
			t2=`cat $file | grep "System time" | awk '{print $4}'`
			tt=`echo "$t1 $t2" | awk '{print $1 + $2}'`
			cc="$cc$tt "
		done
		echo $id $scallop $aa $cc
	done

	# stringtie
	for aa in `echo "tophat star hisat"`
	do
		cc=""
		for abd in `echo "0.01 1 2.5 5 7.5 10 25 50 75 100"`
		do
			file="../$id.$aa/$stringtie.$abd/time.log"

			t1=`cat $file | grep "user" | awk '{print $2}'`
			x1=`echo $t1 | cut -f 1 -d "m"`
			y1=`echo $t1 | cut -f 2 -d "m" | cut -f 1 -d "s"`

			t2=`cat $file | grep "sys" | awk '{print $2}'`
			x2=`echo $t2 | cut -f 1 -d "m"`
			y2=`echo $t2 | cut -f 2 -d "m" | cut -f 1 -d "s"`

			tt=`echo "$x1 $x2 $y1 $y2" | awk '{print $1 * 60 + $2 * 60 + $3 + $4}'`
			
#echo $id stringtie $aa $abd $x1 $x2 $y1 $y2 $tt
			cc="$cc$tt "
		done
		echo $id $stringtie $aa $cc
	done

	# transcomb
	for aa in `echo "tophat star"`
	do
		cc=""
		for abd in `echo "0.01 1 2.5 5 7.5 10 25 50 75 100"`
		do
			file="../$id.$aa/$transcomb.$abd/time.log"

			t1=`cat $file | grep "user" | awk '{print $2}'`
			x1=`echo $t1 | cut -f 1 -d "m"`
			y1=`echo $t1 | cut -f 2 -d "m" | cut -f 1 -d "s"`

			t2=`cat $file | grep "sys" | awk '{print $2}'`
			x2=`echo $t2 | cut -f 1 -d "m"`
			y2=`echo $t2 | cut -f 2 -d "m" | cut -f 1 -d "s"`

			tt=`echo "$x1 $x2 $y1 $y2" | awk '{print $1 * 60 + $2 * 60 + $3 + $4}'`
			cc="$cc$tt "
		done
		echo $id $transcomb $aa $cc
	done

done
