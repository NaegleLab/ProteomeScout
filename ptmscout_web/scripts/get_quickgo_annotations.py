from DB_init import DatabaseInitialization
from ptmscout.database import DBSession, protein
from paste.deploy.loadwsgi import appconfig
import time, datetime
from ptmworker import upload_helpers, quickgo_tools, data_import
import sys, os
from ptmscout.utils.decorators import rate_limit
import traceback

def create_missing_GO_annotations(GO_annotation_map, protein_ids):
    created_go_entries = []
    created, assigned = 0, 0
    
    complete_go_terms = set()
    for acc in GO_annotation_map:
        for goId in GO_annotation_map[acc]:
            complete_go_terms.add(goId)

    protein_accessions = {}
    uniq_ids = set()
    for acc in protein_ids:
        pid = protein_ids[acc]
        uniq_ids.add(pid)

        paccs = protein_accessions.get(pid, set())
        paccs.add(acc)
        protein_accessions[pid] = paccs

    missing_terms = set()

    for protein_id in uniq_ids:
        go_annotations = {}
        
        for acc in protein_accessions[protein_id]:
            for goId in GO_annotation_map[acc]:
                dateAdded = GO_annotation_map[acc][goId]
                if goId not in go_annotations or go_annotations[goId] < dateAdded:
                    go_annotations[goId] = dateAdded
        
        for goId in go_annotations:
            dateAdded = go_annotations[goId]
            entry = upload_helpers.get_go_annotation(goId, protein_id, dateAdded, complete_go_terms, missing_terms)

            if entry:
                created_go_entries.append(entry)
                created+=1

            assigned+=1
    
    missing_terms = list(missing_terms)
    upload_helpers.query_missing_GO_terms(missing_terms)

    print "Assigned %d terms, created %d terms" % ( assigned, created )
    return created_go_entries

def create_hierarchy_for_missing_GO_annotations(created_entries):
    links = 0
    for entry in created_entries:
        go_term = protein.getGoAnnotationById(entry.goId)
        
        for parent_goId in entry.is_a:
            parent_term = protein.getGoAnnotationById(parent_goId)
            
            if parent_term == None:
                print "Parent term '%s' not found in GO annotations!" % (parent_goId)
            else:
                parent_term.children.append(go_term)
                parent_term.save()
                links += 1
    
    print "Created: %d GO entries with %d edges" % (len(created_entries), links)        


@rate_limit(rate=10)
def get_quickgo_annotations(dbconn, protein_map):
    protein_accessions = set([])
    protein_ids = {}
    for pid in protein_map:
        paccs, _ = protein_map[pid]
        for acc in paccs:
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
    created_entries = create_missing_GO_annotations(go_terms, protein_ids)
    print "Created %d terms" % (len(created_entries))
    create_hierarchy_for_missing_GO_annotations(created_entries)

FLUSH_EVERY = 100

if __name__ == "__main__":
    try:
        settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'production.ini'))
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        i = 0
        print "Getting protein data..."

        if len(sys.argv) > 1:
            query_ids = [int(pid) for pid in sys.argv[1:]]
            proteins = DBSession.query(protein.Protein).filter(protein.Protein.id.in_(query_ids)).all()
        else:
            proteins = DBSession.query(protein.Protein).all()

        protein_map = {}

        queries = []
        for prot in proteins:
            paccessions = [ acc.value for acc in prot.accessions ]
            if paccessions != []:
                goIds = prot.getGOIds()
                protein_map[prot.id] = (paccessions, goIds)
            
            i+=1
            if i % FLUSH_EVERY == 0:
                if len(protein_map) > 0:
                    queries.append(protein_map)
                protein_map = {}

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
