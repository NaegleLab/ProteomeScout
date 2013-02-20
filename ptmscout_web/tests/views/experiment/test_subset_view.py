from pyramid.testing import DummyRequest
from ptmscout.views.experiment import subset_view
from mock import patch
from ptmscout.config import strings
from tests.views.mocking import createMockExperiment
from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase

class ExperimentSubsetIntegrationTest(IntegrationTestCase):
    def test_integration(self):
        self.ptmscoutapp.get('/experiments/26/subsets', status=200)

class ExperimentSubsetViewTests(UnitTestCase):

    @patch('ptmscout.views.dataset.dataset_explorer_view.format_explorer_view')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_subset_view(self, patch_getExperiment, patch_format_explorer):
        mock_experiment = createMockExperiment()
        mock_experiment.measurements = ["some","list","of","measurements"]

        patch_getExperiment.return_value = mock_experiment
        patch_format_explorer.return_value = {"some":"map","of":"data"}
        
        mock_experiment.name = "experiment name"
        
        request = DummyRequest()
        request.user = None
        request.matchdict['id'] = 1 
        
        parameters = subset_view.view_experiment_subset(request)

        patch_format_explorer.assert_called_once_with(mock_experiment.measurements)
        patch_getExperiment.assert_called_once_with(1, request.user)
        self.assertEqual(mock_experiment, parameters['experiment'])
        self.assertEqual(strings.experiment_subset_page_title, parameters['pageTitle'])
        self.assertEqual("map", parameters['some'])
        self.assertEqual("data", parameters['of'])

