#!/bin/bash

num=18
gene="ABI3BP"


prefix="$gene".gr
output="$gene".pdf

rm -rf $output

ss=""
for i in `seq 0 $num`
do
	s="$prefix"."$i".pdf
	ss="$ss $s"
	echo $s
	myepstool.sh "$prefix"."$i" 1> /dev/null 2>/dev/null
done

pdftk $ss cat output $output
rm "$prefix".*


prefix="$gene".nt
output="$gene".x.pdf

rm -rf $output

ss=""
for i in `seq 0 $num`
do
	s="$prefix"."$i".pdf
	ss="$ss $s"
	echo $s
	myepstool.sh "$prefix"."$i" 1> /dev/null 2>/dev/null
done

pdftk $ss cat output $output
rm "$prefix".*
