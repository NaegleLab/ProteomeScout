from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
import json
from ptmscout.views.dataset import dataset_explorer_view
from mock import patch, call
from tests.views.mocking import createMockExperiment, createMockMeasurement,\
    createMockUser, createMockPeptide,\
    createMockPeptideModification, createMockPTM, createMockScansite,\
    createMockDataItem
from ptmscout.database import protein

class TestDatasetExplorerView(UnitTestCase):
    
    def test_parse_query_expression(self):
        foreground_query = [
               ['nop', [ 'protein', 'P03367' ]],
               ['and', [ 'metadata', 'Modification Type', 'eq', 'Phosphoserine']],
               ['or',  [ 'quantitative', '2.0', 'gt', 'run1:data:5', 'div', 'run1:data:1' ]], 
               ['and', [ 'scansite', 'Scansite-Kinase', 'eq', 'PDGFRB_Y_kin', '0.2' ]]
            ]
        
        selector = dataset_explorer_view.parse_query_expression(foreground_query)

        d1 = createMockDataItem('run1', 1, 'data', '0', 'time(sec)', 0.4)
        d2 = createMockDataItem('run1', 2, 'data', '1', 'time(sec)', 0.1)
        d3 = createMockDataItem('run1', 3, 'data', '5', 'time(sec)', 0.01)
        d4 = createMockDataItem('run1', 4, 'stddev', '1', 'time(sec)', 0.001)
        d5 = createMockDataItem('run1', 5, 'stddev', '5', 'time(sec)', 0.002)
        passing_series = [d1,d2,d3,d4,d5]
        
        d1 = createMockDataItem('run1', 1, 'data', '0', 'time(sec)', 0.1)
        d2 = createMockDataItem('run1', 2, 'data', '1', 'time(sec)', 0.2)
        d3 = createMockDataItem('run1', 3, 'data', '5', 'time(sec)', 0.5)
        d4 = createMockDataItem('run1', 4, 'stddev', '1', 'time(sec)', 0.01)
        d5 = createMockDataItem('run1', 5, 'stddev', '5', 'time(sec)', 0.002)
        failing_series = [d1,d2,d3,d4,d5]
        
        d1 = createMockDataItem('run1', 1, 'data', '0', 'time(sec)', 0.1)
        d2 = createMockDataItem('run1', 2, 'data', '1', 'time(sec)', 0.2)
        d3 = createMockDataItem('run1', 3, 'data', '5', 'time(sec)', None)
        d4 = createMockDataItem('run1', 4, 'stddev', '1', 'time(sec)', 0.01)
        d5 = createMockDataItem('run1', 5, 'stddev', '5', 'time(sec)', 0.002)
        missing_values_series = [d1,d2,d3,d4,d5]
         
        
        ptm1 = createMockPTM(name='phosphoserine')
        ptm2 = createMockPTM(name='phosphothreonine')
        ptm3 = createMockPTM(name='phosphotyrosine')
        
        p1 = protein.Protein()
        acc1 = protein.ProteinAccession()
        acc1.type = 'uniprot'
        acc1.value = 'P03367'
        p1.accessions = [acc1]
        p2 = protein.Protein()
        p3 = protein.Protein()
        
        ms1 = createMockMeasurement(p1.id, 101)
        ms1.protein = p1
        createMockPeptideModification(ms1, createMockPeptide(p1.id, '1', 'Y'), ptm3)
        ms1.data = failing_series
        
        ms2 = createMockMeasurement(p1.id, 101)
        ms2.protein = p1
        createMockPeptideModification(ms2, createMockPeptide(p1.id, '31', 'S'), ptm1)
        ms2.data = failing_series
        
        ms3 = createMockMeasurement(p1.id, 101)
        ms3.protein = p1
        
        pep = createMockPeptide(p1.id, '101', 'T')
        scansite = createMockScansite(pep.id)
        scansite.value = 'PDGFRB_Y_kin'
        scansite.source = 'scansite_bind'
        scansite.score = 0.01
        pep.predictions = [scansite]
        createMockPeptideModification(ms3, pep, ptm2)
        ms3.data = passing_series
        
        ms4 = createMockMeasurement(p2.id, 101)
        ms4.protein = p2
        createMockPeptideModification(ms4, createMockPeptide(p2.id, '11', 'S'), ptm1)
        ms4.data = passing_series
        
        ms5 = createMockMeasurement(p2.id, 101)
        ms5.protein = p2
        
        pep = createMockPeptide(p2.id, '1011', 'Y')
        scansite = createMockScansite(pep.id)
        scansite.value = 'PDGFRB_Y_kin'
        scansite.source = 'scansite_kinase'
        scansite.score = 0.01
        pep.predictions = [scansite]
        createMockPeptideModification(ms5, pep, ptm3)
        ms5.data = passing_series
        
        
        ms6 = createMockMeasurement(p2.id, 101)
        ms6.protein = p2
        pep = createMockPeptide(p2.id, '21', 'S')
        scansite = createMockScansite(pep.id)
        scansite.value = 'PDGFRB_Y_kin'
        scansite.source = 'scansite_kinase'
        scansite.score = 0.01
        pep.predictions = [scansite]
        createMockPeptideModification(ms6, pep, ptm1)
        ms6.data = failing_series
        
        ms7 = createMockMeasurement(p3.id, 101)
        ms7.protein = p3
        createMockPeptideModification(ms7, createMockPeptide(p3.id, '43', 'Y'), ptm3)
        ms7.data = failing_series
        
        ms8 = createMockMeasurement(p3.id, 101)
        ms8.protein = p3
        pep = createMockPeptide(p3.id, '101', 'T')
        scansite = createMockScansite(pep.id)
        scansite.value = 'PDGFRB_Y_kin'
        scansite.source = 'scansite_kinase'
        scansite.score = 1.0
        pep.predictions = [scansite]
        createMockPeptideModification(ms8, pep, ptm2)
        ms8.data = passing_series
        
        
        ms9 = createMockMeasurement(p2.id, 101)
        ms9.protein = p2
        pep = createMockPeptide(p2.id, '31', 'S')
        scansite = createMockScansite(pep.id)
        scansite.value = 'PDGFRB_Y_kin'
        scansite.source = 'scansite_kinase'
        scansite.score = 0.01
        pep.predictions = [scansite]
        createMockPeptideModification(ms9, pep, ptm1)
        ms9.data = missing_values_series
        
        result = [ ms for ms in [ms1, ms2, ms3, ms4, ms5, ms6, ms7, ms8] if selector(ms) ]
        
        self.assertEqual([ms2, ms5], result)
        
    
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
        self.assertEqual("a seqlogo", result['foreground']['seqlogo'])
        self.assertEqual(3, result['foreground']['proteins'])
        self.assertEqual(3, result['foreground']['peptides'])
        self.assertEqual("a seqlogo", result['background']['seqlogo'])
        self.assertEqual(3, result['background']['proteins'])
        self.assertEqual(3, result['background']['peptides'])
        self.assertEqual("Subset 1", result['name'])
        
        
class IntegrationTestDatasetExplorerView(IntegrationTestCase):
    def test_view_integration(self):
        exp_id = 28
        foreground_query = [
               ['nop', [ 'metadata', 'GO-Biological Process', 'eq', 'GO:0000187: activation of MAPK activity' ]]
            ]
        
        query_expression = {
                            'experiment': exp_id,
                            'type': 'create',
                            'name': 'Subset 1',
                            'background': 'experiment',
                            'foreground': foreground_query}
        
        result = self.ptmscoutapp.post_json("/webservice/subsets", params=query_expression, status=200)
        result = result.json
        
        self.assertEqual(foreground_query, result['foreground']['query'])
        self.assertEqual(4, result['foreground']['proteins'])
        self.assertEqual(7, result['foreground']['peptides'])
        
        self.assertEqual('experiment', result['background']['query'])
        self.assertEqual(52, result['background']['proteins'])
        self.assertEqual(68, result['background']['peptides'])
        
        
        self.assertEqual('Subset 1', result['name'])
        self.assertEqual(exp_id, result['experiment'])
        