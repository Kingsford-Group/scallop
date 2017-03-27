#!/bin/bash

if [ "$#" != 1 ]
then
	echo "usage: file"
	exit
fi

tmp=./tmp.xxxx.cc

echo "/*" > $tmp
echo "Part of Scallop Transcript Assembler" >> $tmp
echo "(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University." >> $tmp
echo "See LICENSE for licensing." >> $tmp
echo "*/" >> $tmp
echo "" >> $tmp

cat $1 >> $tmp

mv $tmp $1
