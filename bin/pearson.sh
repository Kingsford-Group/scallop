#!/bin/bash

salmon="salmon1"
tag="B504"
sra=~/data/transcriptomics/SRA
bin=`pwd`/

#for k in `echo "SRR534307"`
for k in `cat list10`
do
	gtf=$sra/$k/$salmon.gtf
	for m in `echo "hisat"`
	do
		echo $k $m

		dir=$bin/$k.$m
		cur=$dir/scallop.$tag

		$bin/gtfformat RPKM2TPM $cur/scallop.gtf $cur/scallop.tpm.gtf
		$bin/gtfcuff quant $cur/cuffcmp.scallop.gtf.tmap $cur/scallop.tpm.gtf $gtf > $cur/pearson

		cur=$dir/stringtie.`stringtie --version`
		$bin/gtfcuff quant $cur/cuffcmp.st.gtf.tmap $cur/st0.gtf $gtf > $cur/pearson
	done

	for m in `echo "tophat star"`
	do
		echo $k $m

		dir=$bin/$k.$m

		cur=$dir/scallop.$tag
		$bin/gtfformat RPKM2TPM $cur/scallop.gtf $cur/scallop.tpm.gtf
		$bin/gtfcuff quant $cur/cuffcmp.scallop.gtf.tmap $cur/scallop.tpm.gtf $gtf > $cur/pearson

		cur=$dir/stringtie.`stringtie --version`
		$bin/gtfcuff quant $cur/cuffcmp.st.gtf.tmap $cur/st0.gtf $gtf > $cur/pearson

		cur=$dir/transcomb
		$bin/gtfformat FPKM2TPM $cur/TransComb0.gtf $cur/TransComb.tpm.gtf
		$bin/gtfcuff quant $cur/cuffcmp.TransComb.gtf.tmap $cur/TransComb.tpm.gtf $gtf > $cur/pearson

	done
done
