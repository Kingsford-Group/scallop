#!/bin/bash

../scallop ../config expression.sim.gtf > log

ss=""
for i in `cat hard.list | cut -f 1 -d " "`
do
	echo $i
	myepstool.sh $i.gr.0 1> /dev/null 2> /dev/null
	myepstool.sh $i.nt.0 1> /dev/null 2> /dev/null
	ss="$ss $i.gr.0.pdf $i.nt.0.pdf"
done
pdftk $ss cat output example.hard.pdf
rm -rf *.0.*
