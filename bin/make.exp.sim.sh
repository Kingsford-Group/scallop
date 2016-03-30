#!/bin/bash

num=10
dir=./simulation1
params=flux.params
gtf=GSM981256/genes1.gtf

for i in `seq 1 $num`
do
	echo "simulate $i"
#flux-simulator -p $params -x
	cur=$dir/$i
#mkdir -p $cur
#mv profile $cur
	./merge.expression.pl $gtf $cur/profile > $cur/expression.gtf
done
