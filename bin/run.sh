#!/bin/bash

tag="B505.A"
abd=$1

list="GSM981256.star,first \
	  GSM981256.hisat,first \
	  GSM981244.star,first \
	  GSM981244.hisat,first \
	  GSM984609.star,first \
	  GSM984609.hisat,first \
	  SRR307911.star,first \
	  SRR307911.hisat,first \
	  SRR387661.star,first \
	  SRR387661.hisat,first \
	  GSM981256.tophat,first \
	  GSM981244.tophat,first \
	  GSM984609.tophat,first \
	  SRR307911.tophat,first \
	  SRR387661.tophat,first \
	  SRR307903.hisat,first \
	  SRR307903.star,first \
	  SRR307903.tophat,first \
	  SRR315323.hisat,first \
	  SRR315323.star,first \
	  SRR315323.tophat,first \
	  SRR315334.hisat,first \
	  SRR315334.star,first \
	  SRR315334.tophat,first \
	  SRR534307.hisat,first \
	  SRR534307.star,first \
	  SRR534307.tophat,first \
	  SRR545723.hisat,first \
	  SRR545723.star,first \
	  SRR545723.tophat,first"

for i in `echo $list`
do
	id=`echo $i | cut -f 1 -d ","`
	strand=`echo $i |cut -f 2 -d ","`
	aligner=`echo $id | cut -f 2 -d "."`
	echo $id $strand $aligner

#nohup ./run.stringtie.sh $id $aligner "1.3.2d.A.$abd" "-c $abd --rf" &
#nohup ./run.scallop.sh $id $aligner "$tag.$abd" "--min_transcript_coverage $abd --library_type $strand" &
#nohup ./run.scallop.sh $id $aligner "$tag" "--library_type $strand" &
#nohup ./run.stringtie.sh $id $aligner `stringtie --version` &
done

list="GSM981256.star,first \
	  GSM981244.star,first \
	  GSM984609.star,first \
	  SRR307911.star,first \
	  SRR387661.star,first \
	  GSM981256.tophat,first \
	  GSM981244.tophat,first \
	  GSM984609.tophat,first \
	  SRR307911.tophat,first \
	  SRR387661.tophat,first \
	  SRR307903.star,first \
	  SRR307903.tophat,first \
	  SRR315323.star,first \
	  SRR315323.tophat,first \
	  SRR315334.star,first \
	  SRR315334.tophat,first \
	  SRR534307.star,first \
	  SRR534307.tophat,first \
	  SRR545723.star,first \
	  SRR545723.tophat,first"


for i in `echo $list`
do
	id=`echo $i | cut -f 1 -d ","`
	strand=`echo $i |cut -f 2 -d ","`
	aligner=`echo $id | cut -f 2 -d "."`

	nohup ./run.transcomb.sh $id $aligner "A.$abd" "-s $strand -f $abd" &
#nohup ./run.transcomb.sh $id $aligner &
done
