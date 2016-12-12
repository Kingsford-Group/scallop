#!/bin/bash

if [ "$#" != 1 ]
then
	exit
fi

gene=$1

prefix="$gene"
output="$gene".pdf

curdir=`pwd`
texdir=$curdir/tex

rm -rf $output

ss=""
for i in `seq 0 1000`
do
	if [ ! -f $texdir/$prefix.$i.tex ]
	then
		break
	fi
	s="$prefix"."$i".pdf
	ss="$ss $texdir/$s"
	echo $s
	cd $texdir
	myepstool.sh "$prefix"."$i" 1> /dev/null 2>/dev/null
	cd $curdir
done

pdftk $ss cat output $output
#rm "$prefix".*

exit

prefix="$gene"
output="$gene".nt.pdf

rm -rf $output

ss=""
for i in `seq 0 1000`
do
	if [ ! -f $texdir/$prefix.$i.nt.tex ]
	then
		break
	fi

	s="$prefix"."$i".nt.pdf
	ss="$ss $texdir/$s"
	echo $s
	cd $texdir
	myepstool.sh "$prefix"."$i".nt 1> /dev/null 2>/dev/null
	cd $curdir
done

pdftk $ss cat output $output
#rm "$prefix".*
