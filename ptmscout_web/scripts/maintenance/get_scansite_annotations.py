from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, modifications, protein
from ptmworker.helpers import upload_helpers, scansite_tools
from ptmscout.utils import decorators
from ptmworker import peptide_tasks
from collections import defaultdict
import os, sys
import traceback
import datetime

FLUSH_EVERY = 10

def get_motif_class(taxonomy):
    motif_class = None
    if 'mammalia' in taxonomy:
        motif_class="MAMMALIAN"
    elif 'saccharomycotina' in taxonomy:
        motif_class="YEAST"
    elif 'saccharomyces' in taxonomy:
        motif_class="YEAST"

    return motif_class

@decorators.page_query()
def page_proteins(limit, offset):
    sq = DBSession.query(protein.Protein.id).order_by(protein.Protein.id).limit(limit).offset(offset).subquery()
    return DBSession.query(protein.Protein).join(sq, protein.Protein.id == sq.c.id)

def iterate( iterable, process_func, report_func, report_interval=1000 ):
    i = 0
    for obj in iterable:
        process_func(obj)
        i += 1
        if i % report_interval == 0:
            report_func(i)


def assign_scansite_to_protein(prot, scansite_predictions):
    k = 0
    for sp in scansite_predictions:
        pred = protein.ProteinScansite()
        pred.score = sp.score
        pred.percentile = sp.percentile
        pred.value = sp.nickname
        pred.source = sp.parse_source()
        pred.site_pos = int(sp.site[1:])
        
        if not prot.hasPrediction(pred.source, pred.value, pred.site_pos):
            k += 1
            prot.scansite.append(pred)
        
    prot.saveProtein()
    return k

if __name__ == "__main__":
    try:
        try:
            settings = sys.argv[1]
        except:
            settings = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')


        print "Loading settings:",settings
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        total_proteins = DBSession.query(protein.Protein).count()

        j=0
        k=0
        print "Getting scansite annotations..."

        def process_protein(prot):
            global j, k
            motif_class = get_motif_class( upload_helpers.get_taxonomic_lineage( prot.species.name ) )

            if motif_class:
                scansite_predictions = scansite_tools.get_scansite_protein_motifs(prot.sequence, motif_class)
                j += len(scansite_predictions)
                k += assign_scansite_to_protein(prot, scansite_predictions)

        def report_progress(i):
            global j, k
            print "Got %d annotations, assigned %d annotations, processed %d / %d protein sequences" % (j, k, i, total_proteins)
            DBSession.flush()
            dbinit.commit()
            dbinit.new_transaction()

        iterate( page_proteins, process_protein, report_progress, report_interval=FLUSH_EVERY )

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()

        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
