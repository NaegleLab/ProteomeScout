from pyramid import testing
import unittest
from pyramid.testing import DummyRequest
from ptmscout.experiments import experiment_listing
from mock import patch

class ExperimentViewTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
    @patch('ptmscout.database.experiment.getAllExperiments')
    def test_experiment_listing(self, patch_getAllExperiments):
        request = DummyRequest()
        patch_getAllExperiments.return_value = [1,2,3,4,5]
        
        result = experiment_listing(request)
        
        self.assertEqual([1,2,3,4,5], result['experiments'])
        self.assertEqual("Home - Experiments", result['pageTitle'])
        
        