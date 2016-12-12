#!/bin/bash

dir=./EP.star
id=star_without_gene_models

mouse=./genomes/mouse.mm9.genome
map1=$dir/read.transcript.map
map2=$dir/transcript.gene.map
bamfile=$dir/$id.bam
samfile=$dir/$id.sam
gtffile=$dir/expression.gtf

sam0=$dir/sam0
sam=$dir/sam
bam=$dir/bam
gtf=$dir/gtf

mkdir -p $sam0
mkdir -p $sam
mkdir -p $bam
mkdir -p $gtf

samtools view $bamfile > $samfile

./split.sam.pl $samfile $map1 $map2 $sam0

for i in `ls $sam0`
do
	./add.sam.header.pl $sam0/$i $mouse > $sam/$i	
done

for i in `ls $sam`
do
	j=${i/.sam/}
	samtools view -b $sam/$i > $bam/$j.bam
done

./split.gtf.pl $gtffile $gtf
