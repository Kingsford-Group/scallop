#!/bin/bash

dir=`pwd`/salmon
exp=$dir/expression-data-sparse
out=$dir/expression.data2

for i in `ls $exp`
do
	j=${i/.txt/}
	echo $j
	./map.salmon.exp.pl $dir/all_transcript_names.txt $exp/$i > $out/$j
done
