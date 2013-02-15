from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, experiment
from ptmworker.helpers import upload_helpers
from ptmscout.config import settings
from geeneus import Proteome
import traceback
import sys, os



FLUSH_FREQ = 1000
if __name__ == "__main__":
    try:
        config_options = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'test.ini')
        DatabaseInitialization.setUpClass(config_options)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        loaded_exps = [exp for exp in DBSession.query(experiment.Experiment) if exp.status=='loaded']
        unannotated_exps = [ exp for exp in loaded_exps if len(exp.measurements) > 0 and len(exp.modifications) == 0 ]

        for exp in unannotated_exps:
            print "Annotating %d" % (exp.id)

            residues_modified, final_ptms = upload_helpers.summarize_experiment_load(exp.measurements)

            print "Residues: ", len(residues_modified), " -- ", ','.join( residues_modified )
            print "PTMS:     ", len(final_ptms), " -- ", ','.join( [ ptm.name for ptm in final_ptms ] )

            exp.modified_residues = "".join(residues_modified)
            for ptm in final_ptms:
                exp.modifications.append(ptm)

            DBSession.add(exp)

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
