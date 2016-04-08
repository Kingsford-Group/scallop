#!/bin/bash

if [ "$#" != 1 ]
then
	exit
fi

gene=$1

prefix="$gene".gr
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

prefix="$gene".nt
output="$gene".x.pdf

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
