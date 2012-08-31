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
        
    @patch('ptmscout.database.experiment.getExperimentTree')
    def test_experiment_listing(self, patch_getExperimentTree):
        request = DummyRequest()
        patch_getExperimentTree.return_value = [1,2,3,4,5]
        
        result = experiment_listing(request)
        
        self.assertEqual([1,2,3,4,5], result['experiments'])
        self.assertEqual("Home - Experiments", result['pageTitle'])
        
        