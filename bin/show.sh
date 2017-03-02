#!/bin/bash

sc=./
st=../stringtie/

tmp1=$sc/tmp1
tmp2=$sc/tmp2
tmp3=$sc/tmp3
tmp4=$sc/tmp4

cat $sc/cmp2.log | grep TRUE | cut -f 3 -d " " | sort > $tmp1
cat $st/cmp2.log | grep TRUE | cut -f 3 -d " " | sort > $tmp2

comm -13 $tmp1 $tmp2 | sort > $tmp3
cat $st/cmp2.log | grep TRUE | sort -k3,3 > $tmp4

join $tmp3 $tmp4 -11 -23 | sort -k4,4n -k7,7

rm $tmp1 $tmp2 $tmp3 $tmp4
