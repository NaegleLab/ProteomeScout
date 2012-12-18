from DB_init import DatabaseInitialization
from ptmscout.database import DBSession, protein
from paste.deploy.loadwsgi import appconfig
import time, datetime
from ptmworker import upload_helpers, quickgo_tools, data_import
import sys, os
from ptmscout.utils.decorators import rate_limit
import traceback

@rate_limit(rate=10)
def get_quickgo_annotations(dbconn, protein_map):
    protein_accessions = set([])
    protein_ids = {}
    for pid in protein_map:
        acc, _ = protein_map[pid]
        protein_ids[acc] = pid
        protein_accessions.add(acc)

    total_duplicated = 0
    new_terms = 0
    go_terms = data_import.get_GO_annotations(protein_accessions)
    for acc in go_terms:
        duplicated = set()
        for goId in go_terms[acc]:
            if goId in protein_map[protein_ids[acc]][1]:
                duplicated.add(goId)

        for goId in duplicated:
            del go_terms[acc][goId]
        total_duplicated += len(duplicated)
        new_terms += len(go_terms[acc])

    print "Adding %d term refs, not adding %d duplicates" % (new_terms, total_duplicated)

    created_entries = data_import.create_missing_GO_annotations(go_terms, protein_ids)
    dbconn.new_transaction()

    print "Created %d terms" % (len(created_entries))
    data_import.create_hierarchy_for_missing_GO_annotations(created_entries)
    dbconn.new_transaction()

def get_primary_accession(protein, used):
    for acc in protein.accessions:
        if acc.type == "swissprot" and acc.value not in used:
            return acc.value
    for acc in protein.accessions:
        if acc.type == "refseq" and acc.value not in used:
            return acc.value

FLUSH_EVERY = 100

if __name__ == "__main__":
    try:
        settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'production.ini'))
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        i = 0
        print "Getting protein data..."

        proteins = DBSession.query(protein.Protein).all()

        used = set()
        protein_map = {}

        queries = []
        for prot in proteins:
            pacc = get_primary_accession(prot, used)

            if pacc != None:
                goIds = prot.getGOIds()
                protein_map[prot.id] = (pacc, goIds)
                used.add(pacc)
            
            i+=1
            if i % FLUSH_EVERY == 0:
                print protein_map
                if len(protein_map) > 0:
                    queries.append(protein_map)
                protein_map = {}
                used = set()

                print i, "/", len(proteins)
       
        if len(protein_map) > 0:
            queries.append(protein_map)

        i=1
        print "Performing quickGo queries..."
        for protein_map in queries:
            i+=1
            print i, "/", len(queries)
            get_quickgo_annotations(dbinit, protein_map)

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
