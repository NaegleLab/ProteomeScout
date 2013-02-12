import sys, os
from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, uniprot, modifications, experiment
from ptmworker.helpers import upload_helpers
import traceback

FLUSH_EVERY=10

def delete_ambiguities(ms):
    for amb in ms.ambiguities:
        DBSession.delete(amb)
    ms.ambiguities = []

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

        ambiguous_experiments = []
        for exp in DBSession.query(experiment.Experiment):
            if exp.ambiguity == 1 and exp.status == 'loaded':
                ambiguous_experiments.append(exp.id)
                print "Found experiment: %d -- %s" % (exp.id, exp.name)

        i = 0
        j = 0
        print "Processing peptides..."
        for exp_id in ambiguous_experiments:
            measured_peps = modifications.getMeasuredPeptidesByExperiment(exp_id, secure=False, check_ready=False)

            for ms in measured_peps:
                delete_ambiguities(ms)
                upload_helpers.check_ambiguity(ms, ms.protein.species.name)
                j += len(ms.ambiguities)
                DBSession.add(ms)

                i+=1
                if i % FLUSH_EVERY == 0:
                    print "%d processed" % (i)
                    DBSession.flush()

        print "Annotated %d ambiguities" % (j)

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
