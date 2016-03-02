#!/bin/bash

if [ "$#" != 1 ]
then
	echo "usage: $0 <gtf-file>"
	exit
fi

tmp1=/tmp/gtf.tmp1
tmp2=/tmp/gtf.tmp2

#cat $1 | awk '$3 == "exon"' $1 | sed 's/Curated Genomic/Curated_Genomic/g' > $tmp1
cat $1 | sed 's/Curated Genomic/Curated_Genomic/g' > $tmp1

awk 'BEGIN {FS="\t"; OFS="\t"}
{
	split($NF,a," ");
	pfx="";
	s="";
	for(i=1;i<=length(a);i+=2)
	{
		if(a[i]=="transcript_id") {pfx=a[i]" "a[i+1]}
		else {s=s" "a[i]" "a[i+1]}
	}
	if(pfx=="") {print "[WARN] line "NR" without transcript_id!" > "/dev/stderr"}
	else {$NF=pfx""s;print$0} 
}' $tmp1 
