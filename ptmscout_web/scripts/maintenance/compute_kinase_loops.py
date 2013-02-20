import sys, os
from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, protein
from ptmworker.helpers import upload_helpers
from ptmscout.config import strings
import traceback

FLUSH_EVERY = 100

def delete_kinase_regions(prot):
    removable = [ strings.kinase_loop_name, strings.possible_kinase_name ]
    remove = []
    for r in prot.regions:
        if r.label in removable:
            remove.append(r)

    for r in remove:
        prot.regions.remove(r)
        DBSession.delete(r)
    
    return len(remove) > 0

if __name__ == "__main__":
    database = sys.argv[1]

    if database == 'test':
        dbconfig = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'test.ini')
    elif database == 'production':
        dbconfig = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')


    try:
        DatabaseInitialization.setUpClass(dbconfig)
        dbinit = DatabaseInitialization()
        dbinit.setUp()
        
        j = 0
        i = 0
        prot_cnt = DBSession.query(protein.Protein.id).count()

        print "Computing kinase loops..."
        for prot in DBSession.query(protein.Protein):
            deleted = delete_kinase_regions(prot)
            j+=upload_helpers.find_activation_loops(prot)

            if j > 0 or deleted:
                DBSession.add(prot)

            i+=1
            if i % FLUSH_EVERY == 0:
                DBSession.flush()
                print "%d / %d" % (i, prot_cnt)

        print "Annotated %d kinase loops" % (j)
        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
