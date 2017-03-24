#!/bin/bash

if [ "$#" != "1" ]
then
	echo "run.sh abd"
	exit
fi

abd=$1

dir=/home/mingfus/data/repositories/scallop/bin/ENCODE/
list=$dir/alignments.unique.rainforest.list

cmd=commands
rm -rf $cmd

IFS="
"

tag1="B668"
tag2="1.3.2d"
for x in `cat $list`
do
	id=`echo $x | cut -f 1 -d " "`
	gm=`echo $x | cut -f 5 -d " "`
	echo $id $gm

	file1=$dir/scallop.$tag1.$abd/$id/gffmul.stats
	file2=$dir/stringtie.$tag2.$abd/$id/gffmul.stats

#echo "./run.scallop.encode.sh $id $gm $tag1.$abd '--min_transcript_coverage $abd'" >> $cmd
#echo "./run.stringtie.encode.sh $id $gm $tag2.$abd '-c $abd'" >> $cmd

	if [ ! -s $file1 ]; then 
		nohup ./run.scallop.encode.sh $id $gm $tag1.$abd "--min_transcript_coverage $abd" &
	else
		ss=`cat $file1 | grep Intron | grep chain | wc -l`
		if [ "$ss" == "0" ]; then
			nohup ./run.scallop.encode.sh $id $gm $tag1.$abd "--min_transcript_coverage $abd" &
		fi
	fi

	if [ ! -s $file2 ]; then 
		nohup ./run.stringtie.encode.sh $id $gm $tag2.$abd "-c $abd" &
	else
		ss=`cat $file2 | grep Intron | grep chain | wc -l`
		if [ "$ss" == "0" ]; then
			nohup ./run.stringtie.encode.sh $id $gm $tag2.$abd "-c $abd" &
		fi
	fi

#nohup ./run.transcomb.encode.sh $id $gm $abd "-s unstranded -f $abd" &
done
