#!/bin/bash

cd /data/ptmscout/ptmscout_web
python scripts/export/export_compendia.py > data/export/ptmscout_everything.tsv
python scripts/export/export_compendia.py --species "vertebrata" > data/export/ptmscout_vertebrata.tsv
python scripts/export/export_compendia.py --species "mammalia" > data/export/ptmscout_mammalia.tsv
python scripts/export/export_compendia.py --modification "PHOS" > data/export/ptmscout_phosphorylation.tsv
python scripts/export/export_compendia.py --modification "Acetylation" > data/export/ptmscout_acetylation.tsv
python scripts/export/export_compendia.py --modification "METH" > data/export/ptmscout_methylation.tsv
python scripts/export/export_compendia.py --modification "UBIQ" > data/export/ptmscout_ubiquitination.tsv
python scripts/export/export_compendia.py --modification "Glycosylation" > data/export/ptmscout_glycosylation.tsv

python scripts/export/summarize_compendia.py ptmscout_everything.tsv ptmscout_vertebrata.tsv ptmscout_mammalia.tsv ptmscout_phosphorylation.tsv ptmscout_acetylation.tsv ptmscout_methylation.tsv ptmscout_ubiquitination.tsv ptmscout_glycosylation.tsv > data/export/listing.pyp
