from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, modifications
from ptmworker.helpers import upload_helpers, scansite_tools
from ptmworker import peptide_tasks
from collections import defaultdict
import os
import traceback
import datetime

def getMotifClass(taxonomy):
    motif_class = None
    if 'mammalia' in taxonomy:
        motif_class="MAMMALIAN"
    elif 'saccharomycotina' in taxonomy:
        motif_class="YEAST"
    elif 'saccharomyces' in taxonomy:
        motif_class="YEAST"

    return motif_class

FLUSH_EVERY=10

def copyScansitePredictions(pep, scansite_predictions):
    db_predictions = []

    for scansite in scansite_predictions:
        pred = modifications.ScansitePrediction()
        pred.score = scansite.score
        pred.percentile = scansite.percentile
        pred.value = scansite.nickname
        pred.source = scansite.parse_source()
        db_predictions.append(pred)

    pep.predictions = db_predictions
    pep.scansite_date = datetime.datetime.now()


if __name__ == "__main__":
    try:
        settings = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        print "Getting peptide data..."
        unprocessed_peptides = [ pep for pep in DBSession.query(modifications.Peptide).filter(modifications.Peptide.scansite_date==None) ]
        peptide_queries = defaultdict(set)
        sequence_to_peptides = {'MAMMALIAN':defaultdict(set), 'YEAST':defaultdict(set)}
        pep_map = {}

        def map_sum(m):
            i = 0
            for k in m:
                i += len(m[k])
            return i

        i=0
        j=0
        print "Checking for unqueried peptides..."
        for pep in unprocessed_peptides:
            pep_map[pep.id] = pep
            species = pep.protein.species.name
            taxonomy = upload_helpers.get_taxonomic_lineage(species)
            motif_class = getMotifClass(taxonomy)

            if motif_class != None:
                peptide_queries[motif_class].add(pep.pep_aligned)
                sequence_to_peptides[motif_class][pep.pep_aligned].add(pep.id)

            i+=1
            if i % FLUSH_EVERY == 0:
                print "Processed %d / %d found %d unqueried unique peptide sequences" % (i, len(unprocessed_peptides), map_sum(peptide_queries))

        total_peptides = len(peptide_queries['MAMMALIAN']) + len(peptide_queries['YEAST'])
        i = 0
        j = 0
        k = 0
        print "Getting scansite annotations..."
        for motif_class in peptide_queries:
            for pep_seq in peptide_queries[motif_class]:
                scansite_predictions = scansite_tools.get_scansite_motif(pep_seq, motif_class)
                j += len(scansite_predictions)

                for pep_id in sequence_to_peptides[motif_class][pep_seq]:
                    copyScansitePredictions(pep_map[pep_id], scansite_predictions)
                    DBSession.add(pep_map[pep_id])
                    k += 1

                i+=1
                if i % FLUSH_EVERY == 0:
                    print "Got %d annotations, assigned %d / %d peptides, processed %d / %d peptide sequences" % (j, k, len(unprocessed_peptides), i, total_peptides)
                    DBSession.flush()
                    dbinit.commit()
                    dbinit.new_transaction()
            
            print "Loaded %d scansite predictions" % (j)

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
