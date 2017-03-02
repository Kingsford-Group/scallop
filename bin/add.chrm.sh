#!/bin/bash

for i in `cat list10`
do
	echo $i

	dir=$i.hisat/scallop.B504
	mv $dir/scallop.gtf $dir/scallop_no.gtf
	cat $dir/scallop_no.gtf | sed 's/^/chr/g' > $dir/scallop.gtf

	continue

	dir=$i.hisat/stringtie.1.3.1c
	mv $dir/st.gtf $dir/st_no.gtf
	cat $dir/st_no.gtf | sed 's/^/chr/g' > $dir/st.gtf

	dir=$i.hisat/transcomb
	mv $dir/TransComb.gtf $dir/TransComb_no.gtf
	cat $dir/TransComb_no.gtf | sed 's/^/chr/g' > $dir/TransComb.gtf
done
