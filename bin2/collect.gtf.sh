#!/bin/bash

if [ "$#" -ne 2 ];
then
	echo "usage $0 directory full/STRG"
	exit
fi

script=./scritps.sh
rm -rf $script

touch $script
chmod +x $script

for i in `ls $1 | grep gtf`
do
	k=${i/.gtf/}
	echo "cat $1/$i | sed 's/$2\./$k./g'" >> $script
done

./$script
