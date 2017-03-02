#!/bin/bash

tag="B414"

list="CIDANE.tophat \
	  EP.tophat \
	  EP.star \
	  ER.tophat \
	  ER.star \
	  GSM981256.tophat \
	  GSM981256.star \
	  GSM981256.hisat \
	  GSM981244.tophat \
	  GSM981244.star \
	  GSM981244.hisat \
	  GSM984609.tophat \
	  GSM984609.star \
	  GSM984609.hisat \
	  SRR387662.tophat \
	  SRR387662.star \
	  SRR387662.hisat \
	  SRR307911.tophat \
	  SRR307911.star \
	  SRR307911.hisat"

for i in `echo $list`
do
	aligner=`echo $i | cut -f 2 -d "."`
	echo $i $aligner
	nohup ./run.scallop.sh $i $aligner $tag &
	nohup ./run.stringtie.sh $i $aligner `stringtie --version` &
done

list="CIDANE.tophat,second \
	  EP.tophat,second \
	  EP.star,second \
	  ER.tophat,second \
	  ER.star,second \
	  GSM981256.tophat,first \
	  GSM981256.star,first \
	  GSM981244.tophat,first \
	  GSM981244.star,first \
	  GSM984609.tophat,first \
	  GSM984609.star,first \
	  SRR387662.tophat,first \
	  SRR387662.star,first \
	  SRR307911.tophat,first \
	  SRR307911.star,first"

for i in `echo $list`
do
	id=`echo $i | cut -f 1 -d ","`
	strand=`echo $i |cut -f 2 -d ","`
	aligner=`echo $id | cut -f 2 -d "."`
	echo $id $aligner $strand
	nohup ./run.transcomb.sh $id $aligner &
done
