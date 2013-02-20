from ptmscout.database import DBSession
import os, sys
import traceback
from ptmscout.database import protein
from scripts.DB_init import DatabaseInitialization
from obo_parser import GeneOntologyFile

if __name__ == '__main__':
    
    settings = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'development.ini')
    DatabaseInitialization.setUpClass(settings)
    dbinit = DatabaseInitialization()

    try:
        FLUSH_FREQUENCY = 1000
        
        dbinit.setUp()
        
        GO_file_location = "scripts/data/gene_ontology.1_2.obo"
        GO_file = GeneOntologyFile(GO_file_location)

        print "Loading OBO:"
        print "version: ", GO_file.format_version
        print "date:    ", str(GO_file.date)

        term = GO_file.next_term()

        print "inserting terms."
        i = 0
        new_go_terms_created = 0
        while term != None:
            if protein.getGoAnnotationById(term.id) == None:
                dbterm = protein.GeneOntology()
                dbterm.GO = term.id
                dbterm.aspect = term.get_aspect()
                dbterm.term = term.name
                dbterm.date = GO_file.date
                dbterm.version = GO_file.format_version
                dbterm.save()
                new_go_terms_created += 1

            i += 1
            if i % FLUSH_FREQUENCY == 0:
                print ".",
                sys.stdout.flush()
                DBSession.flush()

            term = GO_file.next_term()
        
        print "Created %d GO terms" % (new_go_terms_created)

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
