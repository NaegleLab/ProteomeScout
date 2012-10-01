import unittest
from pyramid import testing
from ptmscout.views.experiment.summary_view import experiment_summary_view
from pyramid.testing import DummyRequest
from ptmscout.config import strings
from tests.views.mocking import createMockExperiment
from mock import patch

class SummaryViewsTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        
    def tearDown(self):
        testing.tearDown()
        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_summary_view(self, patch_getExperiment):
        request = DummyRequest()
        request.matchdict['id'] = 1
        request.user = None
        
        exp = createMockExperiment(1, 0, 0)
        patch_getExperiment.return_value = exp
        
        result = experiment_summary_view(request)
        
        patch_getExperiment.assert_called_once_with(1, None)
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.experiment_summary_page_title % (exp.name), result['pageTitle'])