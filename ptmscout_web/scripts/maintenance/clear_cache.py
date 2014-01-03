from scripts.DB_init import DatabaseInitialization
from ptmscout.config import settings
from ptmscout.database import DBSession, jobs
import sys, os
import traceback

if __name__ == "__main__":
    try:
        config_options = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        DatabaseInitialization.setUpClass(config_options)
        dbinit = DatabaseInitialization()
        dbinit.setUp()
    
        mode = 'clean'
        try:
            mode = sys.argv[1]
        except:
            pass

        cached = DBSession.query(jobs.CachedResult)

        for cr in cached:
            if cr.is_expired() or mode == 'clear':
                print "Delete cache entry (", cr.started, "):", cr.function
                cr.delete()

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
