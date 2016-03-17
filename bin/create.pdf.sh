#!/bin/bash

prefix="ABI3BP.gr"
output="ABI3BP.pdf"

rm -rf $output

for i in `seq 0 18`
do
	myepstool.sh "$prefix"."$i"
	s="$prefix"."$i".pdf
	ss="$ss $s"
done

pdftk $ss cat output $output
rm "$prefix".*
