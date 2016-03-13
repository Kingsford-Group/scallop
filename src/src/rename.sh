#!/bin/bash

tmp=/tmp/rename.tmp
for i in `ls *.h *.cc`
do
	cat $i | sed 's/descriptor_b/descriptor/g' > $tmp; mv $tmp $i
	cat $i | sed 's/iterator_b/iterator/g' > $tmp; mv $tmp $i
	cat $i | sed 's/PEE_b/PEE/g' > $tmp; mv $tmp $i
	cat $i | sed 's/PEB_b/PEB/g' > $tmp; mv $tmp $i
done
