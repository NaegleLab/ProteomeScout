from scripts.DB_init import DatabaseInitialization
from scripts import progressbar
from ptmscout.database import DBSession, modifications, protein
import os, sys
import traceback

if __name__ == "__main__":
    try:
        settings = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', sys.argv[1])
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        total_peptides = DBSession.query(modifications.Peptide).count()
        i = 0
        pb = progressbar.ProgressBar(max_value = total_peptides)
        pb.start()
        for pep in DBSession.query(modifications.Peptide):
            for pred in pep.predictions:
                ps = protein.ProteinScansite()
                
                ps.percentile = pred.percentile
                ps.value = pred.value
                ps.source = pred.source
                ps.score = pred.score
                ps.site_pos = pep.site_pos
                
                if not pep.protein.hasPrediction(ps.source, ps.value, ps.site_pos):
                    pep.protein.scansite.append(ps)
                
                pep.protein.saveProtein()
        
            i+=1
            pb.update(i)
            
        pb.finish()
        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
