from pyramid import testing
import unittest
from pyramid.testing import DummyRequest
from ptmscout.experiments import experiment_listing, view_experiment
from mock import patch, Mock
from ptmscout import strings
from ptmscout.database.experiment import Experiment
from tests.views.test_user_management import createUserForTest
from pyramid.httpexceptions import HTTPForbidden

class ExperimentViewTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
    @patch('ptmscout.database.experiment.getExperimentTree')
    def test_experiment_listing(self, patch_getExperimentTree):
        request = DummyRequest()
        exp1 = Experiment()
        exp2 = Experiment()
        exp3 = Experiment()
        exp1.id = 1
        exp2.id = 2
        exp3.id = 3
        exp1.public = 1
        exp2.public = 1
        exp3.public = 0
        
        patch_getExperimentTree.return_value = [exp1,exp2]
        request.user = None
        
        result = experiment_listing(request)
        
        self.assertEqual([exp1,exp2], result['experiments'])
        self.assertEqual(strings.experiments_page_title, result['pageTitle'])
    
    @patch('ptmscout.database.experiment.getExperimentTree')
    def test_experiment_listing_should_filter_experiments_without_proper_permissions(self, patch_getExperimentTree):
        request = DummyRequest()
        exp1 = Experiment()
        exp2 = Experiment()
        exp1.id = 1
        exp2.id = 2
        exp1.public = 0
        exp2.public = 0
        
        request.user = createUserForTest("user", "email", "password", 1)
        request.user.experiments = []
        request.user.experiments.append(exp1)
        
        patch_getExperimentTree.return_value = [exp1, exp2]
        
        result = experiment_listing(request)
        
        self.assertEqual([exp1], result['experiments'])
        self.assertEqual(strings.experiments_page_title, result['pageTitle'])
        
        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_view(self, patch_getExperiment):
        mock_experiment = Mock('ptmscout.database.experiment.Experiment')
        patch_getExperiment.return_value = mock_experiment
        
        mock_experiment.name = "experiment name"
        mock_experiment.public = 1
        
        request = DummyRequest()
        request.user = None
        request.matchdict['id'] = 1 
        
        parameters = view_experiment(request)
        
        patch_getExperiment.assert_called_once_with(1)
        self.assertEqual(mock_experiment, parameters['experiment'])
        self.assertEqual(strings.experiment_page_title % ("experiment name"), parameters['pageTitle'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_view_should_allow_user_with_permissions(self, patch_getExperiment):
        mock_experiment = Mock('ptmscout.database.experiment.Experiment')
        patch_getExperiment.return_value = mock_experiment
        mock_experiment.id = 1 
        mock_experiment.name = "experiment name"
        mock_experiment.public = 0
        
        request = DummyRequest()
        request.user = createUserForTest("user", "email", "password", 1)
        request.user.experiments = [mock_experiment]
        request.matchdict['id'] = 1 
        
        parameters = view_experiment(request)
        
        patch_getExperiment.assert_called_once_with(1)
        self.assertEqual(mock_experiment, parameters['experiment'])
        self.assertEqual(strings.experiment_page_title % ("experiment name"), parameters['pageTitle'])
    
    
    @patch('ptmscout.database.experiment.getExperimentById')    
    def test_experiment_view_should_raise_forbidden(self, patch_getExperiment):
        mock_experiment = Mock('ptmscout.database.experiment.Experiment')
        patch_getExperiment.return_value = mock_experiment 
        mock_experiment.name = "experiment name"
        mock_experiment.public = 0
        
        request = DummyRequest()
        request.user = None
        request.matchdict['id'] = 1 
        
        try:
            view_experiment(request)
        except HTTPForbidden:
            pass
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected HTTPForbidden")

    @patch('ptmscout.database.experiment.getExperimentById')    
    def test_experiment_view_should_raise_forbidden_for_logged_in_user(self, patch_getExperiment):
        mock_experiment = Mock('ptmscout.database.experiment.Experiment')
        patch_getExperiment.return_value = mock_experiment 
        
        mock_experiment.id = 1
        mock_experiment.name = "experiment name"
        mock_experiment.public = 0
        
        request = DummyRequest()
        request.user = createUserForTest("user", "email", "password", 1)
        request.user.experiments = []
        request.matchdict['id'] = 1 
        
        try:
            view_experiment(request)
        except HTTPForbidden:
            pass
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected HTTPForbidden")
        