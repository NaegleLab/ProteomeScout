from tests.DBTestCase import DBTestCase
import ptmscout.database.experiment as dbexperiment
from ptmscout.database.experiment import NoSuchExperiment, Experiment,\
    getExperimentTree, ExperimentAccessForbidden, ExperimentNotAvailable
from mock import patch
from ptmscout.database import user, permissions

class ExperimentTestCase(DBTestCase):
    def test_getExperimentById_should_succeed_on_existing_experiment(self):
        experiment = dbexperiment.getExperimentById(1, None)
        
        self.assertEqual("PhosphoSite: A bioinformatics resource dedicated to physiological protein phosphorylation.", experiment.name)
        self.assertEqual("http://www.phosphosite.org/", experiment.URL)
    
    def test_getExperimentById_should_throw_exception_when_experiment_is_currently_uploading(self):
        try:
            exp = dbexperiment.getExperimentById(1, None)
            exp.public = 0
            exp.ready = 0
            exp.saveExperiment()
            
            ptmuser = user.getUserById(1)
            
            perm = permissions.Permission(exp)
            perm.experiment_id = exp.id
            perm.user_id = ptmuser.id
            
            ptmuser.permissions.append(perm)
            
            dbexperiment.getExperimentById(1, ptmuser)
            
        except ExperimentNotAvailable, f:
            self.assertEqual(1, f.eid)
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception ExperimentNotAvailable was not thrown")
    
    def test_getExperimentById_should_throw_exception_when_experiment_is_private(self):
        try:
            exp = dbexperiment.getExperimentById(1, None)
            exp.public = 0
            exp.saveExperiment()
            dbexperiment.getExperimentById(1, None)
            
        except ExperimentAccessForbidden, f:
            self.assertEqual(1, f.eid)
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception NoSuchExperiment was not thrown")
        
    
    def test_getExperimentById_should_throw_exception_when_no_such_experiment(self):
        try:
            dbexperiment.getExperimentById(100000, None)
        except NoSuchExperiment, nse:
            self.assertEqual(100000, nse.eid)
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception NoSuchExperiment was not thrown")
            
    def test_getAllExperiments_should_return_all_experiments_in_database_to_which_access_is_allowed(self):
        exp = dbexperiment.getExperimentById(1, None)
        exp.public = 0
        exp.saveExperiment()
        experiments = dbexperiment.getAllExperiments(None)
        
        self.assertTrue(len(experiments) > 0)
        
        self.assertFalse( exp.id in [ e.id for e in experiments ] )
        
    @patch('ptmscout.database.experiment.getAllExperiments')
    def test_getExperimentTree_should_build_experiment_tree(self, patch_getAllExperiments):
        experiments = [Experiment(), Experiment(), Experiment(), Experiment(), Experiment(), Experiment()]
        
        experiments[0].id = 1
        experiments[1].id = 2
        experiments[2].id = 3
        experiments[3].id = 4
        experiments[4].id = 5
        experiments[5].id = 6
        
        experiments[0].experiment_id = None
        experiments[1].experiment_id = None
        experiments[2].experiment_id = 1
        experiments[3].experiment_id = 1
        experiments[4].experiment_id = None
        experiments[5].experiment_id = 2
        
        patch_getAllExperiments.return_value = experiments
        
        current_user = 0
        experiment_tree = getExperimentTree(current_user)
        
        patch_getAllExperiments.assert_called_with(current_user)
        self.assertEqual([experiments[0], experiments[1], experiments[4]], experiment_tree)
        self.assertEqual([experiments[2], experiments[3]], experiments[0].children)
        self.assertEqual([experiments[5]], experiments[1].children)
        