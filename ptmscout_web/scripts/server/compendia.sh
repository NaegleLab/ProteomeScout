#!/bin/bash
set -e

source /data/pyramid/bin/activate

cd /data/ptmscout/ptmscout_web
python scripts/export/export_compendia.py data/export/proteomescout_everything
python scripts/export/export_compendia.py --species "vertebrata" data/export/proteomescout_vertebrata
python scripts/export/export_compendia.py --species "mammalia" data/export/proteomescout_mammalia
python scripts/export/export_compendia.py --modification "PHOS" data/export/proteomescout_phosphorylation
python scripts/export/export_compendia.py --modification "Acetylation" data/export/proteomescout_acetylation
python scripts/export/export_compendia.py --modification "METH" data/export/proteomescout_methylation
python scripts/export/export_compendia.py --modification "UBIQ" data/export/proteomescout_ubiquitination
python scripts/export/export_compendia.py --modification "Glycosylation" data/export/proteomescout_glycosylation

python scripts/export/summarize_compendia.py proteomescout_everything.zip proteomescout_vertebrata.zip proteomescout_mammalia.zip proteomescout_phosphorylation.zip proteomescout_acetylation.zip proteomescout_methylation.zip proteomescout_ubiquitination.zip proteomescout_glycosylation.zip > data/export/listing.pyp

python scripts/maintenance/generate_statistics.py
