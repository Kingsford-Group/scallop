#!/bin/bash

if [ "$#" -ne 3 ] && [ "$#" -ne 4 ]
then
	echo "usage $0 dataset aligner version [parameters]"
	exit
fi

dir=`pwd`/$1/scallop.$3

mkdir -p $dir

bam=`pwd`/$1/$2/"$2".sort.bam
gtf=`pwd`/$1/expression.gtf

{ /usr/bin/time -v ./scallop -i $bam -o $dir/scallop.gtf $4 > $dir/scallop.log; } 2> $dir/time.log
#{ time ./scallop -i $bam -o $dir/scallop.gtf $4 > /dev/null; } 2> $dir/time.log

if [ "$2" = "hisat" ]; then
	mv $dir/scallop.gtf $dir/scallop_no.gtf
	cat $dir/scallop_no.gtf | sed 's/^/chr/g' > $dir/scallop.gtf
fi

cd $dir
gffcompare -r $gtf $dir/scallop.gtf -o gffall
gffcompare -r $gtf $dir/scallop.gtf -o gffmul -M -N
cd -

refallsize=`cat $dir/gffall.stats | grep Reference | grep mRNA | awk '{print $5}'`
refmulsize=`cat $dir/gffmul.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g' `
./gtfcuff roc $dir/gffall.scallop.gtf.tmap $refallsize > $dir/gffall.roc
./gtfcuff roc $dir/gffmul.scallop.gtf.tmap $refmulsize > $dir/gffmul.roc

#./gtfcuff classify $dir/gffmul.scallop.gtf.tmap $dir/scallop.gtf > $dir/gffmul.class

#refmulsize=`cat $dir/gffmul.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g' `
#./gtfcuff roc-trunc $dir/gffmul.scallop.gtf.tmap $refmulsize 0 10.0 > $dir/gffmul.trunc.10

#sra=/home/mingfus/data/transcriptomics/SRA
#id=`echo $1 | cut -f 1 -d "."`
#quant=$sra/"$id".all/salmon/salmon.quant/quant.sf
#./gtfcuff acc-quant $dir/gffmul.scallop.gtf.tmap $quant 0.1 > $dir/gffmul.quant.acc

#refmulsize=`cat $dir/gffmul.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g' `
#./gtfcuff acc $dir/gffmul.scallop.gtf.tmap $refmulsize > $dir/gffmul.acc


exit

sra=~/data/transcriptomics/SRA
genome=/home/mingfus/data/transcriptomics/iGenome/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
k=`echo $1 | cut -f 1 -d "."`

gffread -w $dir/scallop.fa -g $genome $dir/scallop.gtf
salmon index -t $dir/scallop.fa -i $dir/salmon.index 1> /dev/null 2> /dev/null
salmon quant -i $dir/salmon.index -l ISR -1 $sra/$k/$k.Rd1.fq -2 $sra/$k/$k.Rd2.fq -p 4 -o $dir/salmon.quant 1> /dev/null 2> /dev/null
