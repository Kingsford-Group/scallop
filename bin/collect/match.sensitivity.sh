#!/bin/bash

if [ "$#" != "0" ]; then
	echo "usage $0"
	exit
fi

list=../list10
abd="default"

scallop="scallop.B676"
stringtie="stringtie.1.3.2d"
transcomb="transcomb"

for id in `cat $list`
do
	# tophat
	aa="tophat"
	fx="../$id.$aa/$scallop.$abd/gffmul.stats"
	fy="../$id.$aa/$stringtie.$abd/gffmul.stats"
	fz="../$id.$aa/$transcomb.$abd/gffmul.stats"

	tx="../$id.$aa/$scallop.$abd/gffmul.scallop.gtf.tmap"
	ty="../$id.$aa/$stringtie.$abd/gffmul.st.gtf.tmap"
	tz="../$id.$aa/$transcomb.$abd/gffmul.TransComb.gtf.tmap"

	sx=`cat $fx | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g'`
	sy=`cat $fy | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g'`
	sz=`cat $fz | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g'`

	x1=`cat $fx | grep Matching | grep intron | grep chain | awk '{print $4}'`
	y1=`cat $fy | grep Matching | grep intron | grep chain | awk '{print $4}'`
	z1=`cat $fz | grep Matching | grep intron | grep chain | awk '{print $4}'`

	x2=`cat $fx | grep Intron | grep chain | awk '{print $6}'`
	y2=`cat $fy | grep Intron | grep chain | awk '{print $6}'`
	z2=`cat $fz | grep Intron | grep chain | awk '{print $6}'`

	xx1=$x1
	xx2=$x2
	yy1=$y1
	yy2=$y2
	zz1=$z1
	zz2=$z2

	if (( $(echo "$x1 <= $y1" | bc -l) )) && (( $(echo "$x1 <= $z1" | bc -l) )); then
		yy1=`gtfcuff match-correct $ty $sy $x1 | cut -f 10 -d " "`
		zz1=`gtfcuff match-correct $tz $sz $x1 | cut -f 10 -d " "`
		yy2=`gtfcuff match-correct $ty $sy $x1 | cut -f 16 -d " "`
		zz2=`gtfcuff match-correct $tz $sz $x1 | cut -f 16 -d " "`
	elif (( $(echo "$y1 <= $x1" | bc -l) )) && (( $(echo "$y1 <= $z1" | bc -l) )); then
		xx1=`gtfcuff match-correct $tx $sx $y1 | cut -f 10 -d " "`
		xx1=`gtfcuff match-correct $tx $sx $y1 | cut -f 10 -d " "`
		zz1=`gtfcuff match-correct $tz $sz $y1 | cut -f 10 -d " "`
		xx2=`gtfcuff match-correct $tx $sx $y1 | cut -f 16 -d " "`
		zz2=`gtfcuff match-correct $tz $sz $y1 | cut -f 16 -d " "`
	elif (( $(echo "$z1 <= $x1" | bc -l) )) && (( $(echo "$z1 <= $y1" | bc -l) )); then
		xx1=`gtfcuff match-correct $tx $sx $z1 | cut -f 10 -d " "`
		yy1=`gtfcuff match-correct $ty $sy $z1 | cut -f 10 -d " "`
		xx2=`gtfcuff match-correct $tx $sx $z1 | cut -f 16 -d " "`
		yy2=`gtfcuff match-correct $ty $sy $z1 | cut -f 16 -d " "`
	fi

	cc1="$xx1 $xx2 $yy1 $yy2 $zz1 $zz2"
#echo "$id $aa | $x1 $xx1 $x2 $xx2 | $y1 $yy1 $y2 $yy2 | $z1 $zz1 $z2 $zz2 |"

	# star
	aa="star"
	fx="../$id.$aa/$scallop.$abd/gffmul.stats"
	fy="../$id.$aa/$stringtie.$abd/gffmul.stats"
	fz="../$id.$aa/$transcomb.$abd/gffmul.stats"

	tx="../$id.$aa/$scallop.$abd/gffmul.scallop.gtf.tmap"
	ty="../$id.$aa/$stringtie.$abd/gffmul.st.gtf.tmap"
	tz="../$id.$aa/$transcomb.$abd/gffmul.TransComb.gtf.tmap"


	sx=`cat $fx | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g'`
	sy=`cat $fy | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g'`
	sz=`cat $fz | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g'`

	x1=`cat $fx | grep Matching | grep intron | grep chain | awk '{print $4}'`
	y1=`cat $fy | grep Matching | grep intron | grep chain | awk '{print $4}'`
	z1=`cat $fz | grep Matching | grep intron | grep chain | awk '{print $4}'`

	x2=`cat $fx | grep Intron | grep chain | awk '{print $6}'`
	y2=`cat $fy | grep Intron | grep chain | awk '{print $6}'`
	z2=`cat $fz | grep Intron | grep chain | awk '{print $6}'`

	xx1=$x1
	xx2=$x2
	yy1=$y1
	yy2=$y2
	zz1=$z1
	zz2=$z2

	if (( $(echo "$x1 <= $y1" | bc -l) )) && (( $(echo "$x1 <= $z1" | bc -l) )); then
		yy1=`gtfcuff match-correct $ty $sy $x1 | cut -f 10 -d " "`
		zz1=`gtfcuff match-correct $tz $sz $x1 | cut -f 10 -d " "`
		yy2=`gtfcuff match-correct $ty $sy $x1 | cut -f 16 -d " "`
		zz2=`gtfcuff match-correct $tz $sz $x1 | cut -f 16 -d " "`
	elif (( $(echo "$y1 <= $x1" | bc -l) )) && (( $(echo "$y1 <= $z1" | bc -l) )); then
		xx1=`gtfcuff match-correct $tx $sx $y1 | cut -f 10 -d " "`
		zz1=`gtfcuff match-correct $tz $sz $y1 | cut -f 10 -d " "`
		xx2=`gtfcuff match-correct $tx $sx $y1 | cut -f 16 -d " "`
		zz2=`gtfcuff match-correct $tz $sz $y1 | cut -f 16 -d " "`
	elif (( $(echo "$z1 <= $x1" | bc -l) )) && (( $(echo "$z1 <= $y1" | bc -l) )); then
		xx1=`gtfcuff match-correct $tx $sx $z1 | cut -f 10 -d " "`
		yy1=`gtfcuff match-correct $ty $sy $z1 | cut -f 10 -d " "`
		xx2=`gtfcuff match-correct $tx $sx $z1 | cut -f 16 -d " "`
		yy2=`gtfcuff match-correct $ty $sy $z1 | cut -f 16 -d " "`
	fi

#echo $id $aa $x1 $xx1 $x2 $xx2 $y1 $yy1 $y2 $yy2 $z1 $zz1 $z2 $zz2
	cc2="$xx1 $xx2 $yy1 $yy2 $zz1 $zz2"


	# hisat
	aa="hisat"
	fx="../$id.$aa/$scallop.$abd/gffmul.stats"
	fy="../$id.$aa/$stringtie.$abd/gffmul.stats"

	tx="../$id.$aa/$scallop.$abd/gffmul.scallop.gtf.tmap"
	ty="../$id.$aa/$stringtie.$abd/gffmul.st.gtf.tmap"

	sx=`cat $fx | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g'`
	sy=`cat $fy | grep Reference | grep mRNA | awk '{print $9}' | sed 's/(//g'`

	x1=`cat $fx | grep Matching | grep intron | grep chain | awk '{print $4}'`
	y1=`cat $fy | grep Matching | grep intron | grep chain | awk '{print $4}'`

	x2=`cat $fx | grep Intron | grep chain | awk '{print $6}'`
	y2=`cat $fy | grep Intron | grep chain | awk '{print $6}'`

	xx1=$x1
	xx2=$x2
	yy1=$y1
	yy2=$y2

	if (( $(echo "$x1 <= $y1" | bc -l) )); then
		yy1=`gtfcuff match-correct $ty $sy $x1 | cut -f 10 -d " "`
		yy2=`gtfcuff match-correct $ty $sy $x1 | cut -f 16 -d " "`
	else
		xx1=`gtfcuff match-correct $tx $sx $y1 | cut -f 10 -d " "`
		xx2=`gtfcuff match-correct $tx $sx $y1 | cut -f 16 -d " "`
	fi

	cc3="$xx1 $xx2 $yy1 $yy2"

	echo $id $cc1 $cc2 $cc3
done

