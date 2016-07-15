#!/bin/bash

name=p2
data=`pwd`/../data/ensembl/$1/gtf/
idmap=$data/id.map
gtf=$data/$name.gtf
sgtf=$data/"$name"_sorted.gtf

dir=`pwd`/ensembl/$1/exp1
num=10

mkdir -p $dir

for i in `seq 1 $num`
do
	cur=$dir/$i
	mkdir -p $cur
	ln -sf $gtf $cur/
	ln -sf $sgtf $cur/
	params=$cur/params

	echo $i

#./merge.sim.exp.pl $gtf $cur/params.pro > $cur/expression.gtf
	./split.bed.pl $cur/params.bed $idmap $cur/genes
	continue;

	echo "REF_FILE_NAME	$name.gtf" > $params
	echo "GEN_DIR         /home/mingfus/data/transcriptomics/iGenomes/Homo.sapiens/GRCh38/Sequence/Chromosomes" >> $params
	echo "NB_MOLECULES    5000000" >> $params
	echo "TSS_MEAN	50" >> $params
	echo "POLYA_SCALE     100" >> $params
	echo "POLYA_SHAPE     1.5" >> $params
	echo "FRAG_SUBSTRATE	RNA" >> $params
	echo "FRAG_METHOD	UR" >> $params
	echo "FRAG_UR_ETA     350" >> $params
	echo "RTRANSCRIPTION	YES" >> $params
	echo "RT_MOTIF	default" >> $params
	echo "PCR_DISTRIBUTION default" >> $params
	echo "GC_MEAN 	 NaN" >> $params
	echo "PCR_PROBABILITY  0.05" >> $params
	echo "FILTERING 	NO" >> $params
	echo "READ_NUMBER	150000000" >> $params
	echo "READ_LENGTH	75" >> $params
	echo "PAIRED_END	YES" >> $params
#flux-simulator -p $params -x -l -s
done
