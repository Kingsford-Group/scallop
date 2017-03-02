#!/bin/bash

dir=`pwd`/
for k in `cat list`
do
	echo $k
	cur=$dir/$k
	mkdir -p $cur

	scp -r mingfus@arctic:$cur/scallop.B502.B $cur
	scp -r mingfus@arctic:$cur/stringtie.1.3.1c $cur
	scp -r mingfus@arctic:$cur/transcomb.v.1.0 $cur
	cd $cur
	ln -sf transcomb.v.1.0 transcomb
	cd -
done
