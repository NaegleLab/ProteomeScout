#!/bin/bash

rootfile=$1
i=0
openfile=$rootfile.$i
MAX_LINES=$2

echo "" > $openfile	
header=`head -n1 $rootfile`

j=0
while read line
do
	echo "$line" >> $openfile
	let j=j+1
	if [ $j -eq $MAX_LINES ]
	then
		let i=i+1
		let j=0
		openfile=$rootfile.$i
		echo "$header" > $openfile
	fi
done < $rootfile

