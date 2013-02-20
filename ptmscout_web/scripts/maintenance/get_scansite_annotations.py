from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, modifications
from ptmworker.helpers import upload_helpers
from ptmworker import peptide_tasks
import os
import traceback

FLUSH_EVERY=100

if __name__ == "__main__":
    try:
        settings = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        valid_taxons = set(['mammalia', 'saccharomyces', 'saccharomycotina'])

        print "Getting peptide data..."
        unprocessed_peptides = [ pep for pep in DBSession.query(modifications.Peptide) if pep.scansite_date==None ]
        unqueried_peptides = []

        i=0
        print "Checking for unqueried peptides..."
        for pep in unprocessed_peptides:
            species = pep.protein.species.name
            taxonomy = upload_helpers.get_taxonomic_lineage(species)

            if len(set(taxonomy) & valid_taxons)>0:
                unqueried_peptides.append(pep)

            i+=1
            if i % FLUSH_EVERY == 0:
                print "Processed %d / %d found %d unqueried peptides" % (i, len(unprocessed_peptides), len(unqueried_peptides))

        i = 0
        j = 0
        print "Getting scansite annotations..."
        for pep in unqueried_peptides:
            species = pep.protein.species.name
            taxonomy = upload_helpers.get_taxonomic_lineage(species)
            peptide_tasks.load_scansite_peptide(pep, taxonomy)

            j+=len(pep.predictions)

            i+=1
            if i % FLUSH_EVERY == 0:
                print "Added %d annotations, processed %d / %d peptides " % (j, i, len(unqueried_peptides))
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
