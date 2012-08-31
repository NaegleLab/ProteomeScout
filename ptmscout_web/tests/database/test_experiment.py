from tests.DBTestCase import DBTestCase
import ptmscout.database.experiment as dbexperiment
from ptmscout.database.experiment import NoSuchExperiment, Experiment,\
    getExperimentTree
from mock import patch

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
        
    @patch('ptmscout.database.experiment.getAllExperiments')
    def test_getExperimentTree_should_build_experiment_tree(self, patch_getAllExperiments):
        experiments = [Experiment(), Experiment(), Experiment(), Experiment(), Experiment(), Experiment()]
        
        experiments[0].id = 1
        experiments[1].id = 2
        experiments[2].id = 3
        experiments[3].id = 4
        experiments[4].id = 5
        experiments[5].id = 6
        
        experiments[0].experiment_id = 0
        experiments[1].experiment_id = 0
        experiments[2].experiment_id = 1
        experiments[3].experiment_id = 1
        experiments[4].experiment_id = 0
        experiments[5].experiment_id = 2
        
        patch_getAllExperiments.return_value = experiments
                
        experiment_tree = getExperimentTree()
                
        self.assertEqual([experiments[0], experiments[1], experiments[4]], experiment_tree)
        self.assertEqual([experiments[2], experiments[3]], experiments[0].children)
        self.assertEqual([experiments[5]], experiments[1].children)
        