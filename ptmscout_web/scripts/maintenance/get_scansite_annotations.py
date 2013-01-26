from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, modifications, taxonomies
from paste.deploy.loadwsgi import appconfig
import time, datetime
from ptmworker import upload_helpers
import sys, os
from ptmscout.utils.decorators import rate_limit
import traceback

FLUSH_EVERY=100

if __name__ == "__main__":
    try:
        settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'production.ini'))
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        i = 0
        j = 0
        print "Getting peptide data..."

        peptides = DBSession.query(modifications.Peptide).all()
        for pep in peptides:

            if pep.scansite_date==None:
                species = pep.protein.species.name
                taxonomy = upload_helpers.get_taxonomic_lineage(species)
                peptide_tasks.load_scansite_peptide(pep, taxonomy)
                j+=len(pep.predictions)

            if i % FLUSH_EVERY == 0:
                print i, "/", len(peptides)
                DBSession.flush()
            i+=1
        
        print "Loaded %d scansite predictions" % (j)

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
