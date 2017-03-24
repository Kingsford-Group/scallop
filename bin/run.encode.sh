#!/bin/bash

if [ "$#" != "1" ]
then
	echo "run.sh abd"
	exit
fi

abd=$1

#list=/home/mingfus/data/repositories/scallop/bin/ENCODE/example.list
list=/home/mingfus/data/repositories/scallop/bin/ENCODE/alignments.unique.rainforest.list

IFS="
"

tag="B668"
for x in `cat $list`
do
	id=`echo $x | cut -f 1 -d " "`
	gm=`echo $x | cut -f 5 -d " "`
	echo $id $gm
	nohup ./run.scallop.encode.sh $id $gm $tag.$abd "--min_transcript_coverage $abd" &
done

tag="1.3.2d"
for x in `cat $list`
do
	id=`echo $x | cut -f 1 -d " "`
	gm=`echo $x | cut -f 5 -d " "`
	echo $id $gm
	nohup ./run.stringtie.encode.sh $id $gm $tag.$abd "-c $abd" &
done

for x in `cat $list`
do
	id=`echo $x | cut -f 1 -d " "`
	gm=`echo $x | cut -f 5 -d " "`
	echo $id $gm
	nohup ./run.transcomb.encode.sh $id $gm $abd "-s unstranded -f $abd" &
done
