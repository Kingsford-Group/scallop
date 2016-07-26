#!/bin/bash

tmp0=/tmp/tmp0
tmp1=/tmp/tmp1
tmp2=/tmp/tmp2
tmp3=/tmp/tmp3
tmp4=/tmp/tmp4
tmp5=/tmp/tmp5

ls $1 | grep -v gtf | grep -v log | grep -v track | grep -v loci > $tmp0

s=`cat $tmp0`

cd $1

cat $s | grep mRNA | grep Query | cut -c 26-27 > $tmp1
cat $s | grep mRNA | grep Reference | cut -c 26-27 > $tmp2

paste $tmp0 $tmp1 $tmp2 | awk '$2 != 0' > $tmp3
paste $tmp0 $tmp1 $tmp2 | awk '$2 == 0' > $tmp4

cat $s | grep chain | grep Intron | cut -c 22-26 > $tmp1
cat $s | grep chain | grep Intron | cut -c 28-32 > $tmp2

paste $tmp3 $tmp1 $tmp2
paste $tmp4 $tmp4 $tmp4 | cut -f 1,2,3,5,8
