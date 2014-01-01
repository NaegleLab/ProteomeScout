#!/bin/bash

files=$(ls logs/ptmworker.*.pid 2> /dev/null)

for f in $files
do
    echo "Stopping workers... $f"
    sudo kill `cat "$f"`
    while [ -f "$f" ]
    do
        sleep 1
    done
done

echo "Done."
