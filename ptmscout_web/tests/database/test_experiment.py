from tests.DBTestCase import DBTestCase
import ptmscout.database.experiment as dbexperiment
from ptmscout.database.experiment import NoSuchExperiment

class ExperimentTestCase(DBTestCase):
    def test_getExperimentById_should_succeed_on_existing_experiment(self):
        experiment = dbexperiment.getExperimentById(1)
        
        self.assertEqual("PhosphoSite: A bioinformatics resource dedicated to physiological protein phosphorylation.", experiment.name)
        self.assertEqual("http://www.phosphosite.org/", experiment.URL)
    
    def test_getExperimentById_should_throw_exception_when_no_such_experiment(self):
        try:
            dbexperiment.getExperimentById(100000)
        except NoSuchExperiment, nse:
            self.assertEqual(100000, nse.eid)
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception NoSuchExperiment was not thrown")
            
    def test_getAllExperiments_should_return_all_experiments_in_database(self):
        experiments = dbexperiment.getAllExperiments()
        self.assertTrue(len(experiments) > 0)