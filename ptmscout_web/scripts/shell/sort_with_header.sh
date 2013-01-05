#!/bin/bash

file=$1
head -n1 $file > $file.sorted.txt
tail +2 $file | sort -k2 >> $file.sorted.txt
