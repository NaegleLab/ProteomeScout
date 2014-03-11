#!/bin/bash
set -e

source /data/pyramid/bin/activate

cd /data/ptmscout/ptmscout_web
python scripts/export/export_compendia.py > data/export/proteomescout_everything.tsv
python scripts/export/export_compendia.py --species "vertebrata" > data/export/proteomescout_vertebrata.tsv
python scripts/export/export_compendia.py --species "mammalia" > data/export/proteomescout_mammalia.tsv
python scripts/export/export_compendia.py --modification "PHOS" > data/export/proteomescout_phosphorylation.tsv
python scripts/export/export_compendia.py --modification "Acetylation" > data/export/proteomescout_acetylation.tsv
python scripts/export/export_compendia.py --modification "METH" > data/export/proteomescout_methylation.tsv
python scripts/export/export_compendia.py --modification "UBIQ" > data/export/proteomescout_ubiquitination.tsv
python scripts/export/export_compendia.py --modification "Glycosylation" > data/export/proteomescout_glycosylation.tsv

python scripts/export/summarize_compendia.py proteomescout_everything.tsv proteomescout_vertebrata.tsv proteomescout_mammalia.tsv proteomescout_phosphorylation.tsv proteomescout_acetylation.tsv proteomescout_methylation.tsv proteomescout_ubiquitination.tsv proteomescout_glycosylation.tsv > data/export/listing.pyp

python scripts/maintenance/generate_statistics.py
