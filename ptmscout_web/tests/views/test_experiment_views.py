from pyramid import testing
import unittest
from pyramid.testing import DummyRequest
from ptmscout.experiments import experiment_listing, view_experiment
from mock import patch, Mock
from ptmscout import strings

class ExperimentViewTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
    @patch('ptmscout.database.experiment.getExperimentTree')
    def test_experiment_listing(self, patch_getExperimentTree):
        request = DummyRequest()
        patch_getExperimentTree.return_value = [1,2,3,4,5]
        
        result = experiment_listing(request)
        
        self.assertEqual([1,2,3,4,5], result['experiments'])
        self.assertEqual(strings.experiments_page_title, result['pageTitle'])
        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_view(self, patch_getExperiment):
        mock_experiment = Mock('ptmscout.database.experiment.Experiment')
        patch_getExperiment.return_value = mock_experiment 
        mock_experiment.name = "experiment name"
        
        request = DummyRequest()
        request.matchdict['id'] = 1 
        
        parameters = view_experiment(request)
        
        patch_getExperiment.assert_called_once_with(1)
        self.assertEqual(mock_experiment, parameters['experiment'])
        self.assertEqual(strings.experiment_page_title % ("experiment name"), parameters['pageTitle'])