#!/bin/bash

lines=${2:-40}
tail -n$lines -f $1 | cut -c1-400
