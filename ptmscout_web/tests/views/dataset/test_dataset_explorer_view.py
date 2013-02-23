from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
import json
from ptmscout.views.dataset import dataset_explorer_view
from mock import patch, call
from tests.views.mocking import createMockExperiment, createMockMeasurement,\
    createMockUser

class TestDatasetExplorerView(UnitTestCase):
    
    @patch('ptmscout.views.dataset.dataset_explorer_view.calculate_feature_enrichment')
    @patch('ptmscout.utils.protein_utils.create_sequence_profile')
    @patch('ptmscout.views.dataset.dataset_explorer_view.parse_query_expression')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_subset_query_create_subset(self, patch_getExp, patch_parse, patch_seqlogo, patch_calc):
        exp = createMockExperiment()
        measurements = [createMockMeasurement(2, exp.id),createMockMeasurement(3, exp.id),createMockMeasurement(4, exp.id)]
        exp.measurements = measurements
        
        exp_id = exp.id
        foreground_query = [
                       ['nop', [ 'protein', 'P03367' ]], 
                       ['or',  [ 'quantitative', '2.0', 'gt', 'run1:data:5', 'div', 'run1:data:1' ]], 
                       ['and', [ 'scansite', 'Scansite-Kinase', 'eq', 'PDGFRB_Y_kin', '5' ]]
                    ]
        
        query_expression = {
                            'experiment': exp_id,
                            'type': 'create',
                            'name': 'Subset 1',
                            'background': 'experiment',
                            'foreground': foreground_query}
        
        request = DummyRequest()
        request.user = createMockUser()
        request.body = json.dumps(query_expression)
        
        def comparison_func(ms):
            return True
        
        patch_getExp.return_value = exp
        patch_parse.return_value = comparison_func
        patch_seqlogo.return_value = "a seqlogo"
        patch_calc.return_value = "Some feature enrichment"
        
        result = dataset_explorer_view.compute_subset_enrichment(request)
        
        patch_getExp.assert_called_once_with(exp_id, request.user)
        patch_calc.assert_called_once_with(measurements, measurements)
        
        self.assertIn(call('experiment'), patch_parse.call_args_list)
        self.assertIn(call(foreground_query), patch_parse.call_args_list)
        
        self.assertEqual("Some feature enrichment", result['enrichment'])
        self.assertEqual("a seqlogo", result['foreground'])
        self.assertEqual("a seqlogo", result['background'])
        self.assertEqual("Subset 1", result['name'])
        
        
#class IntegrationTestDatasetExplorerView(IntegrationTestCase):
#    def test_view_integration(self):
#        self.ptmscoutapp.get("/webservice/subsets", status=200)