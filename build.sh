#!/bin/bash

dir=`pwd`

# compile lib
libs="util graph gtf"
for i in `echo $libs`
do
	cur="$dir/lib/$i"
	cd $cur
	aclocal
	autoconf
	autoheader
	automake -a
	./configure
	make
done

# compile src
cd $dir/src
aclocal
autoconf
autoheader
automake -a
./configure
make
