#!/bin/bash

tag="B502.B"
dir=`pwd`/
sra=~/data/transcriptomics/SRA

for k in `cat list10`
do
	for m in `echo "hisat"`
	do
		echo $k $m

		dir=`pwd`/$k.$m

		cur=$dir/scallop.$tag/salmon
		mkdir -p $cur


		cd $cur
		{ time salmon quant -i $index -l ISR -1 $dir/$k.Rd1.fq -2 $dir/$k.Rd2.fq -p 8 -o salmon.quant; } 2> time.log
		cd -

#cur=$dir/scallop.$tag
#cd $cur
#cuffcompare -r $gtf $cur/scallop.gtf
#cd -
#refallsize=`cat $cur/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
#./cuffroc $cur/cuffcmp.scallop.gtf.tmap $cur/scallop.gtf $refallsize > $cur/salmon.all.roc

		cur=$dir/stringtie.`stringtie --version`
		cd $cur
		cuffcompare -r $gtf $cur/st.gtf
		cd -
		refallsize=`cat $cur/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
		./cuffroc $cur/cuffcmp.st.gtf.tmap $cur/st.gtf $refallsize > $cur/salmon.all.roc
	done

	continue

	for m in `echo "tophat star"`
	do
		echo $k $m

		dir=`pwd`/$k.$m

		cur=$dir/scallop.$tag
		cd $cur
		cuffcompare -r $gtf $cur/scallop.gtf
		cd -
		refallsize=`cat $cur/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
		./cuffroc $cur/cuffcmp.scallop.gtf.tmap $cur/scallop.gtf $refallsize > $cur/salmon.all.roc

		cur=$dir/stringtie.`stringtie --version`
		cd $cur
		cuffcompare -r $gtf $cur/st.gtf
		cd -
		refallsize=`cat $cur/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
		./cuffroc $cur/cuffcmp.st.gtf.tmap $cur/st.gtf $refallsize > $cur/salmon.all.roc

		cur=$dir/transcomb
		cd $cur
		cuffcompare -r $gtf $cur/TransComb.gtf
		cd -
		refallsize=`cat $cur/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
		./cuffroc $cur/cuffcmp.TransComb.gtf.tmap $cur/TransComb.gtf $refallsize > $cur/salmon.all.roc
	done

done
