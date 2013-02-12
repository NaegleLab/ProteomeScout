from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, protein
import time, datetime
from ptmworker.helpers import quickgo_tools
from ptmworker import GO_tasks
import sys, os
import traceback

def get_quickgo_annotations(protein_accessions):
    query_accessions = set([ acc for pid in protein_accessions for acc in protein_accessions[pid] ])
    quickgo_result, _ = quickgo_tools.batch_get_GO_annotations(list(query_accessions))

    return quickgo_result

def perform_query_and_insert(protein_accessions, protein_map):
    quickgo_result = get_quickgo_annotations(protein_accessions)
    created_entries = GO_tasks.assign_annotations(quickgo_result, protein_accessions, protein_map)
    GO_tasks.query_missing_GO_terms(created_entries)

    print "Created a total of %d new GO terms" % (len(created_entries))

    GO_tasks.create_hierarchies(created_entries)


def build_queries(proteins):
    protein_accessions = {}
    protein_map = {}

    queries = []
    i = 0
    for prot in proteins:
        paccessions = GO_tasks.get_primary_accessions(prot)

        if paccessions != []:
            protein_accessions[prot.id] = paccessions
            protein_map[prot.id] = prot

        i+=1
        if i % FLUSH_EVERY == 0:
            if len(protein_map) > 0:
                queries.append((protein_accessions, protein_map))
            protein_accessions = {}
            protein_map = {}

            print i, "/", len(proteins)

    if len(protein_map) > 0:
        queries.append((protein_accessions, protein_map))

    return queries


FLUSH_EVERY = 100

if __name__ == "__main__":
    try:
        settings = os.path.join('data', 'ptmscout', 'ptmscout_web', 'production.ini')
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()
        dbinit.log_capture('ptmscout')

        print "Getting protein data..."
        if len(sys.argv) > 1:
            query_ids = [int(pid) for pid in sys.argv[1:]]
            proteins = DBSession.query(protein.Protein).filter(protein.Protein.id.in_(query_ids)).all()
        else:
            proteins = DBSession.query(protein.Protein).all()

        print "Creating queries..."
        queries = build_queries(proteins)

        i=0
        print "Performing quickGo queries..."
        for protein_accessions, protein_map in queries:
            i+=1
            print i, "/", len(queries)
            perform_query_and_insert(protein_accessions, protein_map)

            dbinit.commit()
            dbinit.new_transaction()

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
