from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, modifications, protein
from ptmworker.helpers import upload_helpers, scansite_tools
from ptmscout.utils import decorators
from ptmworker import peptide_tasks
from collections import defaultdict
import os
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

def copy_predictions(pep, scansite_predictions):
    db_predictions = []

    for scansite in scansite_predictions:
        if pep.pep_aligned == scansite.sequence:
            pred = modifications.ScansitePrediction()
            pred.score = scansite.score
            pred.percentile = scansite.percentile
            pred.value = scansite.nickname
            pred.source = scansite.parse_source()
            db_predictions.append(pred)

    pep.predictions = db_predictions
    pep.scansite_date = datetime.datetime.now()

def create_missing_peptides(prot, pep_map, scansite_predictions):
    created = 0
    for scansite in scansite_predictions:
        site_pos = int(scansite.site[1:])
        site_type = scansite.site[0]
        if not any( k == site_pos for k in pep_map ):
            created += 1
            npep = modifications.Peptide()
            npep.site_pos = site_pos
            npep.site_type = site_type
            npep.pep_aligned = prot.get_kmer(site_pos, k=7)
            npep.protein_domain = prot.get_domain(site_pos)
            npep.protein_id = prot.id
            pep_map[site_pos] = npep
    return created

@decorators.page_query()
def page_proteins(limit, offset):
    sq = DBSession.query(protein.Protein.id).order_by(protein.Protein.id).limit(limit).offset(offset).subquery()
    return DBSession.query(protein.Protein).join(sq, protein.Protein.id == sq.c.id)

def get_peptides_for_protein(protein_id):
    pep_map = {}
    for pep in DBSession.query(modifications.Peptide).filter(modifications.Peptide.protein_id==protein_id):
        pep_map[pep.site_pos] = pep
    return pep_map

def assign_scansite_to_peptides(pep_map, scansite_predictions):
    k = 0
    for site in pep_map:
        pep = pep_map[site]
        copy_predictions(pep, scansite_predictions)
        DBSession.add(pep)
        k += 1
    return k

def is_annotated(pep_map):
    return any( len(pep_map[site_pos].predictions) > 0 for site_pos in pep_map )


def iterate( iterable, process_func, report_func, report_interval=1000 ):
    i = 0
    for obj in iterable:
        process_func(obj)
        i += 1
        if i % report_interval == 0:
            report_func(i)

if __name__ == "__main__":
    try:
        settings = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        total_proteins = DBSession.query(protein.Protein).count()

        j=0
        k=0
        created = 0
        print "Getting scansite annotations..."

        def process_protein(prot):
            global j, k, created
            motif_class = get_motif_class( upload_helpers.get_taxonomic_lineage( prot.species.name ) )

            if motif_class:
                pep_map = get_peptides_for_protein(prot.id)
                if not is_annotated(pep_map):
                    scansite_predictions = scansite_tools.get_scansite_protein_motifs(prot.sequence, motif_class)
                    j += len(scansite_predictions)
                    created += create_missing_peptides(prot, pep_map, scansite_predictions)
                    k += assign_scansite_to_peptides(pep_map, scansite_predictions)

        def report_progress(i):
            print "Got %d annotations, assigned %d peptides, created %d peptides, processed %d / %d protein sequences" % (j, k, created, i, total_proteins)
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
