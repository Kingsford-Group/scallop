#!/bin/bash

dir=./ensembl/mouse/

#./purify.gtf.pl $dir/p0.gtf > $dir/p1.gtf
#./fix.gtf.sh $dir/p1.gtf > $dir/p2.gtf

params=flux.exp.params
flux-simulator -p $params -x
mv profile $dir

gtf=$dir/p2.gtf
./merge.exp.pl $gtf $dir/profile > $dir/expression.gtf

./scallop -c ./config -i $dir/expression.gtf -o $dir/s1.gtf -a core -s 5 > $dir/log1
./scallop -c ./config -i $dir/expression.gtf -o $dir/s2.gtf -a full -s 5 > $dir/log2
./scallop -c ./config -i $dir/expression.gtf -o $dir/s3.gtf -a greedy -s 5 > $dir/log3
