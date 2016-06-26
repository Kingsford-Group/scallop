#!/bin/bash

dir=./ensembl/zebrafish2

#./purify.gtf.pl $dir/p0.gtf > $dir/p1.gtf
#./fix.gtf.sh $dir/p1.gtf > $dir/p2.gtf

param0=flux.params.template
params=$dir/params

cp $param0 $params
echo "REF_FILE_NAME	p2.gtf" >> $params
echo "NB_MOLECULES	5000000" >> $params
echo "EXPRESSION_K	-0.1" >> $params

flux-simulator -p $params -x
#mv profile $dir

gtf=$dir/p2.gtf
./merge.exp.pl $gtf $dir/profile > $dir/expression.gtf

./scallop -c ./config -i $dir/expression.gtf -o $dir/s1.gtf -a core -s 5 > $dir/log1
./scallop -c ./config -i $dir/expression.gtf -o $dir/s2.gtf -a full -s 5 > $dir/log2
./scallop -c ./config -i $dir/expression.gtf -o $dir/s3.gtf -a greedy -s 5 > $dir/log3

./gtfcompare $dir/s1.gtf $dir/expression.gtf > $dir/s1.cmp
./gtfcompare $dir/s2.gtf $dir/expression.gtf > $dir/s2.cmp
./gtfcompare $dir/s3.gtf $dir/expression.gtf > $dir/s3.cmp

./gtfformat $dir/expression.gtf $dir/expression2.gtf
