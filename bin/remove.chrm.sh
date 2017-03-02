#!/bin/bash

for i in `cat list10`
do
	echo $i

	dir=$i.hisat/scallop.B502.C
	mv $dir/scallop.gtf $dir/scallop_no.gtf
	cat $dir/scallop_no.gtf | sed 's/^chrchr/chr/g' > $dir/scallop.gtf

	dir=$i.hisat/stringtie.1.3.1c
	mv $dir/st.gtf $dir/st_no.gtf
	cat $dir/st_no.gtf | sed 's/^chrchr/chr/g' > $dir/st.gtf
done
