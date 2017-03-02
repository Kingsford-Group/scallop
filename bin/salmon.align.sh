#!/bin/bash

salmon="salmon"
tag="B504"
dir=`pwd`/
sra=~/data/transcriptomics/SRA
genome=/home/mingfus/data/transcriptomics/iGenome/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa


#for k in `cat list10`
for k in `echo "SRR534307"`
do
	gtf=$sra/$k/$salmon.gtf
	for m in `echo "hisat"`
	do
		echo $k $m

		dir=`pwd`/$k.$m
		cur=$dir/scallop.$tag
		gffread -w $cur/scallop.fa -g $genome $cur/scallop.gtf
		salmon index -t $cur/scallop.fa -i $cur/salmon.index 2> /dev/null
		salmon quant -i $cur/salmon.index -l ISR -1 $sra/$k/$k.Rd1.fq -2 $sra/$k/$k.Rd2.fq -p 8 -o $cur/salmon.quant 2> /dev/null

		cur=$dir/stringtie.`stringtie --version`
	done

	continue;

	for m in `echo "tophat star"`
	do
		echo $k $m

		dir=`pwd`/$k.$m

		cur=$dir/scallop.$tag
		cd $cur
		cuffcompare -r $gtf $cur/scallop.gtf
		cd -
		refallsize=`cat $cur/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
		./cuffroc $cur/cuffcmp.scallop.gtf.tmap $cur/scallop.gtf $refallsize > $cur/$salmon.all.roc

		cur=$dir/stringtie.`stringtie --version`
		cd $cur
		cuffcompare -r $gtf $cur/st.gtf
		cd -
		refallsize=`cat $cur/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
		./cuffroc $cur/cuffcmp.st.gtf.tmap $cur/st.gtf $refallsize > $cur/$salmon.all.roc

		cur=$dir/transcomb
		cd $cur
		cuffcompare -r $gtf $cur/TransComb.gtf
		cd -
		refallsize=`cat $cur/cuffcmp.stats | grep Reference | grep mRNA | awk '{print $5}'`
		./cuffroc $cur/cuffcmp.TransComb.gtf.tmap $cur/TransComb.gtf $refallsize > $cur/$salmon.all.roc
	done

done
