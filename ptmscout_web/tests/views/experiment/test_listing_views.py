from pyramid import testing
import unittest
from pyramid.testing import DummyRequest
from ptmscout.views.experiment.listings_view import experiment_listing, view_experiment
from mock import patch, Mock
from ptmscout.config import strings
from tests.views.mocking import createMockExperiment

class ExperimentListingViewTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
    
    @patch('ptmscout.database.experiment.getExperimentTree')
    def test_experiment_listing(self, patch_getExperimentTree):
        request = DummyRequest()
        exp1 = createMockExperiment(1,1)
        exp2 = createMockExperiment(2,1)
        
        patch_getExperimentTree.return_value = [exp1,exp2]
        request.user = None
        
        result = experiment_listing(request)
        
        patch_getExperimentTree.assert_called_once_with(request.user)
        self.assertEqual([exp1,exp2], result['experiments'])
        self.assertEqual(strings.experiments_page_title, result['pageTitle'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_view(self, patch_getExperiment):
        mock_experiment = Mock('ptmscout.database.experiment.Experiment')
        patch_getExperiment.return_value = mock_experiment
        
        mock_experiment.name = "experiment name"
        
        request = DummyRequest()
        request.user = None
        request.matchdict['id'] = 1 
        
        parameters = view_experiment(request)
        
        patch_getExperiment.assert_called_once_with(1, request.user)
        self.assertEqual(mock_experiment, parameters['experiment'])
        self.assertEqual(strings.experiment_page_title % ("experiment name"), parameters['pageTitle'])
