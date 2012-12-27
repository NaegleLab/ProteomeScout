#!/bin/bash

if [ -e "logs/ptmworker.pid" ]
then
    echo "Stopping workers..."
    kill `cat "logs/ptmworker.pid"`
    while [ -f "logs/ptmworker.pid" ]
    do
        sleep 1
    done
else
    echo "No workers to stop."
fi

echo "Done."
