#!/bin/bash

if [ "$#" != "1" ]
then
	echo "run.sh abd"
	exit
fi

abd=$1

dir=/home/mingfus/data/repositories/scallop/bin/ENCODE/
list=$dir/good.list
#list=$dir/others.list

cmd=commands
rm -rf $cmd

IFS="
"

tag1="B676"
tag2="1.3.2d"
for x in `cat $list`
do
	id=`echo $x | cut -f 1 -d " "`
	ss=`echo $x | cut -f 2 -d " "`
	gm=`echo $x | cut -f 3 -d " "`
	echo $id $ss $gm

	st=""
	if [ "$ss" == "first" ]; then
		st="--rf"
	elif [ "$ss" == "second" ]; then
		st="--fr"
	fi

#nohup ./run.scallop.encode.sh $id $gm $tag1.$abd "--library_type $ss" &
#nohup ./run.stringtie.encode.sh $id $gm $tag2.$abd "$st" &
	nohup ./run.transcomb.encode.sh $id $gm $abd "-s $ss" &

#nohup ./run.scallop.encode.sh $id $gm $tag1.$abd "--min_transcript_coverage $abd --library_type $ss" &
#nohup ./run.stringtie.encode.sh $id $gm $tag2.$abd "-c $abd $st" &
#nohup ./run.transcomb.encode.sh $id $gm $abd "-f $abd -s $ss" &

done
