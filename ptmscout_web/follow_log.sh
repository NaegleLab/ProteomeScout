#!/bin/bash
tail -n20 -f $1 | cut -c1-300
