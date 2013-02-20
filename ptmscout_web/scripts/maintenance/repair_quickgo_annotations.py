from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, protein
from ptmworker.helpers import quickgo_tools
import sys, os
import traceback
import logging
import pickle

FLUSH_EVERY = 100

def update_hierarchy_links(term_entries):
    changed_terms = {}
    links = 0
    i=0
    for goId in term_entries:
        go_term = protein.getGoAnnotationById(goId)
        entry = term_entries[goId]

        for parentId in entry.is_a:
            parent_term = protein.getGoAnnotationById(parentId)

            if not parent_term.hasChild(goId):
                parent_term.children.append(go_term)
                changed_terms[parentId] = parent_term
                links += 1

        i+=1
        if i % FLUSH_EVERY==0:
            print "%d / %d" % (i, len(term_entries))

    print "Created %d new edges" % (links)

    for goId in changed_terms:
        changed_terms[goId].save()


def get_term(goId, missing_terms):
    version, entry = quickgo_tools.get_GO_term(goId)

    for parentId in entry.is_a:
        if protein.getGoAnnotationById(parentId) == None:
            missing_terms.add(parentId)

    return entry, version

def get_all_parents(missing_terms, term_entries):
    created = 0
    while len(missing_terms) > 0:
        goId = missing_terms.pop()
        if goId in term_entries:
            continue

        entry, version = get_term(goId, missing_terms)
        term_entries[goId] = entry

        go_term = protein.GeneOntology()
        go_term.GO = goId
        go_term.term = entry.goName
        go_term.aspect = entry.goFunction
        go_term.version = version
        go_term.save()

        created += 1

    print "Created %d new GO terms" % ( created )


def get_all_terms():
    term_entries = {}
    missing_terms = set()

    go_term_cnt = DBSession.query(protein.GeneOntology).count()
    i = 0
    for go_term in DBSession.query(protein.GeneOntology):
        goId = go_term.GO
        entry, _v = get_term(goId, missing_terms)
        term_entries[goId] = entry

        i+=1
        if i % FLUSH_EVERY == 0:
            print "%d / %d" % (i, go_term_cnt)

    return term_entries, missing_terms

if __name__ == "__main__":
    try:
        settings = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')

        if len(sys.argv) > 1 and sys.argv[1] == '--force':
            os.remove('repair_go.entries')

        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()
        dbinit.log_capture('ptmscout', level=logging.INFO)

        tmp_cache_filename = 'repair_go.entries'

        if os.path.exists(tmp_cache_filename):
            print "Loading terms from cache..."
            with open(tmp_cache_filename, 'r') as cache_file:
                term_entries, missing_terms = pickle.loads(cache_file.read())

            print "Re-checking missingness..."
            still_missing = set()
            for goId in missing_terms:
                if protein.getGoAnnotationById(goId) == None:
                    still_missing.add(goId)
            missing_terms = still_missing

        else:
            print "Loading terms from quickgo..."
            term_entries, missing_terms = get_all_terms()


            print "Saving terms to cache..."
            with open(tmp_cache_filename, 'w') as cache_file:
                cache_file.write(pickle.dumps( (term_entries, missing_terms) ))

        print "Found %d missing parent terms" % ( len(missing_terms) )
        print "Getting missing parent terms..."
        get_all_parents(missing_terms, term_entries)

        print "Updating is_a relationships..."
        update_hierarchy_links(term_entries)

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
