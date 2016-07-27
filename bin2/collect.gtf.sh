#!/bin/bash

script=./scritps.sh
rm -rf $script

touch $script
chmod +x $script

for i in `ls $1 | grep gtf`
do
	k=${i/.gtf/}
	echo "cat $1/$i | sed 's/STRG\./$k./g'" >> $script
done

./$script
