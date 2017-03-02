#!/bin/bash

if [ "$#" -ne 3 ] && [ "$#" -ne 4 ] && [ "$#" -ne 2 ]
then
	echo "usage $0 dataset aligner [version] [parameters]"
	exit
fi

dir=`pwd`/$1/transcomb

if [ "$#" -ge 3 ]
then
	dir=`pwd`/$1/transcomb.$3
fi

mkdir -p $dir

bam=`pwd`/$1/$2/"$2".sort.bam
gtf=`pwd`/$1/expression.gtf

cd $dir
{ /usr/bin/time -v TransComb -b $bam -o $dir $4 > $dir/transcomb.log; } 2> $dir/time.log
cd -

#./gtfcompare $gtf $dir/TransComb.gtf 2 > $dir/cmp.roc

mv $dir/TransComb.gtf $dir/TransComb0.gtf
./gtfformat format $dir/TransComb0.gtf $dir/TransComb.gtf

cd $dir
gffcompare -r $gtf $dir/TransComb.gtf -o gffall
gffcompare -r $gtf $dir/TransComb.gtf -o gffmul -M -N
cd -

refallsize=`cat $dir/gffall.stats | grep Reference | grep mRNA | awk '{print $5}'`
refmulsize=`cat $dir/gffmul.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g'`
./gtfcuff roc $dir/gffall.TransComb.gtf.tmap $refallsize > $dir/gffall.roc
./gtfcuff roc $dir/gffmul.TransComb.gtf.tmap $refmulsize > $dir/gffmul.roc

#./gtfcuff classify $dir/gffmul.TransComb.gtf.tmap $dir/TransComb.gtf > $dir/gffmul.class
#./gtfcuff roc-trunc $dir/gffmul.TransComb.gtf.tmap $refmulsize 0 10.0 > $dir/gffmul.trunc.10

#sra=/home/mingfus/data/transcriptomics/SRA
#id=`echo $1 | cut -f 1 -d "."`
#quant=$sra/"$id".all/salmon/salmon.quant/quant.sf
#./gtfcuff acc-quant $dir/gffmul.TransComb.gtf.tmap $quant 0.1 > $dir/gffmul.quant.acc

#refmulsize=`cat $dir/gffmul.stats | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g'`
#./gtfcuff acc $dir/gffmul.TransComb.gtf.tmap $refmulsize > $dir/gffmul.acc

exit

sra=~/data/transcriptomics/SRA
genome=/home/mingfus/data/transcriptomics/iGenome/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
k=`echo $1 | cut -f 1 -d "."`

gffread -w $dir/TransComb.fa -g $genome $dir/TransComb.gtf
salmon index -t $dir/TransComb.fa -i $dir/salmon.index 1> /dev/null 2> /dev/null
salmon quant -i $dir/salmon.index -l ISR -1 $sra/$k/$k.Rd1.fq -2 $sra/$k/$k.Rd2.fq -p 4 -o $dir/salmon.quant 1> /dev/null 2> /dev/null
