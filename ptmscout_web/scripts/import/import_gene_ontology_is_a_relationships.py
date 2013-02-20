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
        
        missing_terms = set()

        i=0
        term = GO_file.next_term()
        print "inserting relationships."
        while term != None:
            dbterm = protein.getGoAnnotationById(term.id)

            if dbterm == None:
                missing_terms.add(term.id)
            else:
                for pid, pterm in term.is_a_relationships.items():
                    parent_entry = protein.getGoAnnotationById(pid)
                    if parent_entry == None:
                        missing_terms.add(pid)
                    else:
                        parent_entry.children.append(dbterm)
                
            i+=1
            if(i % FLUSH_FREQUENCY == 0):
                print ".",
                sys.stdout.flush()
                DBSession.flush()

            term = GO_file.next_term()
        
        print "Missing terms:", len(missing_terms)
        ofile = open('missing_terms', 'w')
        for term in missing_terms:
            ofile.write(term+"\n")
        ofile.close()
        
        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
