#!/bin/bash
set -e

source /data/pyramid/bin/activate

cd /data/ptmscout/ptmscout_web
python scripts/maintenance/clean_unused_data.py
