#!/bin/bash

cd /data/ptmscout/ptmscout_web
python scripts/export/export_compendia.py > data/export/everything.tsv
python scripts/export/export_compendia.py --species "vertebrata" > data/export/vertebrata.tsv
python scripts/export/export_compendia.py --species "mammalia" > data/export/mammalia.tsv
python scripts/export/export_compendia.py --modification "PHOS" > data/export/phosphorylation.tsv
python scripts/export/export_compendia.py --modification "Acetylation" > data/export/acetylation.tsv
python scripts/export/export_compendia.py --modification "METH" > data/export/methylation.tsv
python scripts/export/export_compendia.py --modification "UBIQ" > data/export/ubiquitination.tsv
python scripts/export/export_compendia.py --modification "Glycosylation" > data/export/glycosylation.tsv

python scripts/export/summarize_compendia.py everything.tsv vertebrata.tsv mammalia.tsv phosphorylation.tsv acetylation.tsv methylation.tsv ubiquitination.tsv glycosylation.tsv > data/export/listing.pyp
