#!/bin/bash

dir=`pwd`

# compile lib
cd $dir/lib
./build.sh

# compile src
cd $dir/src
aclocal
autoconf
autoheader
automake -a
./configure
make
