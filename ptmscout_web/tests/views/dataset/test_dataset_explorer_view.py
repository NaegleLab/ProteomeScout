from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
import json
from ptmscout.views.dataset import dataset_explorer_view
from mock import patch, call, Mock
from tests.views.mocking import createMockExperiment, createMockMeasurement,\
    createMockUser, createMockPeptide,\
    createMockPeptideModification, createMockPTM, createMockScansite,\
    createMockDataItem, createMockSubset
from ptmscout.database import protein, annotations
from ptmscout.config import strings

class TestDatasetExplorerView(UnitTestCase):
    
    def test_parse_query_expression(self):
        foreground_query = [
               ['nop', [ 'protein', 'P03367' ]],
               ['and', [ 'metadata', 'Modification Type', 'eq', 'Phosphoserine']],
               ['or',  [ 'quantitative', '2.0', 'gt', 'run1:data:5', 'div', 'run1:data:1' ]], 
               ['and', [ 'scansite', 'Scansite-Kinase', 'eq', 'PDGFRB_Y_kin', '0.2' ]]
            ]
        
        selector = dataset_explorer_view.parse_query_expression(foreground_query, None, None, {})

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
        scansite = createMockScansite(p1.id, pep.site_pos)
        scansite.value = 'PDGFRB_Y_kin'
        scansite.source = 'scansite_bind'
        scansite.percentile = 0.01
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
        scansite = createMockScansite(p2.id, pep.site_pos)
        scansite.value = 'PDGFRB_Y_kin'
        scansite.source = 'scansite_kinase'
        scansite.percentile = 0.01
        pep.predictions = [scansite]
        createMockPeptideModification(ms5, pep, ptm3)
        ms5.data = passing_series
        
        
        ms6 = createMockMeasurement(p2.id, 101)
        ms6.protein = p2
        pep = createMockPeptide(p2.id, '21', 'S')
        scansite = createMockScansite(p2.id, pep.site_pos)
        scansite.value = 'PDGFRB_Y_kin'
        scansite.source = 'scansite_kinase'
        scansite.percentile = 0.01
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
        scansite = createMockScansite(p2.id, pep.site_pos)
        scansite.value = 'PDGFRB_Y_kin'
        scansite.source = 'scansite_kinase'
        scansite.percentile = 1.0
        pep.predictions = [scansite]
        createMockPeptideModification(ms8, pep, ptm2)
        ms8.data = passing_series
        
        
        ms9 = createMockMeasurement(p2.id, 101)
        ms9.protein = p2
        pep = createMockPeptide(p2.id, '31', 'S')
        scansite = createMockScansite(p2.id, pep.site_pos)
        scansite.value = 'PDGFRB_Y_kin'
        scansite.source = 'scansite_kinase'
        scansite.percentile = 0.01
        pep.predictions = [scansite]
        createMockPeptideModification(ms9, pep, ptm1)
        ms9.data = missing_values_series
        
        result = [ ms for ms in [ms1, ms2, ms3, ms4, ms5, ms6, ms7, ms8] if selector(ms) ]
        
        self.assertEqual([ms2, ms5], result)
        
    
    @patch('ptmscout.views.dataset.dataset_explorer_view.calculate_feature_enrichment')
    @patch('ptmscout.views.dataset.dataset_explorer_view.compute_annotations')
    @patch('ptmscout.utils.protein_utils.create_sequence_profile')
    @patch('ptmscout.views.dataset.dataset_explorer_view.parse_query_expression')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_subset_query_create_subset(self, patch_getExp, patch_parse, patch_seqlogo, patch_annotate, patch_calc):
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
                            'name': 'Subset 1',
                            'background': 'experiment',
                            'foreground': foreground_query,
                            'annotation_set_id': 10}
        
        request = DummyRequest()
        request.user = createMockUser()
        request.body = json.dumps(query_expression)
        request.route_url = Mock()
        request.route_url.return_value = 'http://example.com/proteins/10'
        
        def comparison_func(ms):
            return True
        
        patch_annotate.return_value = {"some": "Annotations"}, {}
        patch_getExp.return_value = exp
        patch_parse.return_value = comparison_func
        patch_seqlogo.return_value = "a seqlogo"
        patch_calc.return_value = "Some feature enrichment"
        
        result = dataset_explorer_view.perform_subset_enrichment(request)
        
        patch_getExp.assert_called_once_with(exp_id, request.user)
        patch_calc.assert_called_once_with(measurements, measurements, {"some": "Annotations"})
        
        self.assertIn(call('experiment', exp_id, request.user, {"some": "Annotations"}), patch_parse.call_args_list)
        self.assertIn(call(foreground_query, exp_id, request.user, {"some": "Annotations"}), patch_parse.call_args_list)
        
        self.assertEqual('success', result['status'])
        self.assertEqual("Some feature enrichment", result['enrichment'])
        self.assertEqual("a seqlogo", result['foreground']['seqlogo'])
        self.assertEqual(3, result['foreground']['proteins'])
        self.assertEqual(3, result['foreground']['peptides'])
        self.assertEqual("a seqlogo", result['background']['seqlogo'])
        self.assertEqual(3, result['background']['proteins'])
        self.assertEqual(3, result['background']['peptides'])
        self.assertEqual("Subset 1", result['name'])
        

    @patch('ptmscout.views.dataset.dataset_explorer_view.compute_annotations')
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.annotations.getSubsetByName')
    def test_save_subset_should_report_error_if_element_exists(self, patch_getSubset, patch_getExp, patch_getAnnotations):
        exp_id = 11
        exp = createMockExperiment(eid=exp_id)

        patch_getSubset.return_value = annotations.Subset()
        patch_getExp.return_value = exp
        
        foreground_query = [
               ['nop', [ 'protein', 'P03367' ]], 
               ['or',  [ 'quantitative', '2.0', 'gt', 'run1:data:5', 'div', 'run1:data:1' ]], 
               ['and', [ 'scansite', 'Scansite-Kinase', 'eq', 'PDGFRB_Y_kin', '5' ]]
            ]
        
        background_query = 'experiment'
        
        save_query = {'name': 'new named subset',
                      'experiment': exp_id,
                      'foreground': foreground_query,
                      'background': background_query}
            
        request = DummyRequest()
        request.user = createMockUser()
        
        request.body = json.dumps(save_query)
        patch_getAnnotations.return_value = {}
        
        result = dataset_explorer_view.save_subset(request)
        
        patch_getExp.assert_called_once_with(exp_id, request.user)
        patch_getSubset.assert_called_once_with(exp_id, 'new named subset', request.user)
        
        self.assertEqual('error', result['status'])
        self.assertEqual(strings.error_message_subset_name_exists, result['message'])
        
    @patch('ptmscout.views.dataset.dataset_explorer_view.compute_annotations')
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.annotations.getSubsetByName')
    @patch('ptmscout.database.annotations.Subset')    
    def test_save_subset_should_create_new_subset_element_and_save(self, patch_subset, patch_getSubset, patch_getExp, patch_getAnnotations):
        exp_id = 11
        exp = createMockExperiment(eid=11)

        patch_getSubset.return_value = None
        patch_getExp.return_value = exp
        
        foreground_query = [
               ['nop', [ 'protein', 'P03367' ]], 
               ['or',  [ 'quantitative', '2.0', 'gt', 'run1:data:5', 'div', 'run1:data:1' ]], 
               ['and', [ 'scansite', 'Scansite-Kinase', 'eq', 'PDGFRB_Y_kin', '5' ]]
            ]
        
        background_query = 'experiment'
        
        save_query = {'name': 'new named subset',
                      'experiment': str(exp_id),
                      'foreground': foreground_query,
                      'background': background_query}
            
        request = DummyRequest()
        request.user = createMockUser()
        user_id = request.user.id
        
        request.body = json.dumps(save_query)
            
        patch_getAnnotations.return_value = {}
        subset = patch_subset.return_value
        
        result = dataset_explorer_view.save_subset(request)
        
        self.assertEqual(exp_id, subset.experiment_id)
        self.assertEqual(user_id, subset.owner_id)
        self.assertEqual('new named subset', subset.name)
        self.assertEqual(foreground_query, subset.foreground_query)
        self.assertEqual(background_query, subset.background_query)
        
        subset.save.assert_called_once_with()
        self.assertEqual('new named subset', result['name'])
        self.assertEqual('success', result['status'])
        
        patch_getExp.assert_called_once_with(exp_id, request.user)
        
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.annotations.getSubsetByName')
    def test_fetch_subset_should_return_error_if_subset_does_not_exist(self, patch_getSubset, patch_getExp):
        exp_id = 11
        exp = createMockExperiment(eid=11)

        patch_getSubset.return_value = None
        patch_getExp.return_value = exp
        
        fetch_query = {
                       'experiment': str(exp_id),
                       'name': 'named subset',
                       'annotation_set_id': 10
                       }
                    
        request = DummyRequest()
        request.user = createMockUser()
        
        request.body = json.dumps(fetch_query)
        
        result = dataset_explorer_view.fetch_subset(request)
        
        self.assertEqual('error', result['status'])
        self.assertEqual(strings.error_message_subset_name_does_not_exist, result['message'])

    @patch('ptmscout.views.dataset.dataset_explorer_view.compute_subset_enrichment')
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.annotations.getSubsetByName')
    def test_fetch_subset_should_return_compute_enrichment_if_subset_exists(self, patch_getSubset, patch_getExp, patch_compute):
        exp_id = 11
        exp = createMockExperiment(eid=11)
        user = createMockUser()
        subset = createMockSubset(name='named subset', owner_id=user.id, experiment_id=exp_id, foreground='some exp', background='some other exp')

        patch_getSubset.return_value = subset
        patch_getExp.return_value = exp
        
        fetch_query = {
                       'experiment': str(exp_id),
                       'name': 'named subset'
                       }
                    
        request = DummyRequest()
        request.user = user
        
        request.body = json.dumps(fetch_query)
        
        patch_compute.return_value = {'some': 'enrichment analysis'}
        
        result = dataset_explorer_view.fetch_subset(request)
        
        self.assertEqual({'id': subset.id, 'some': 'enrichment analysis'}, result)
        patch_compute.assert_called_once_with(request, None, exp, user, 'named subset', 'some exp', 'some other exp', None)
        patch_getExp.assert_called_once_with(exp_id, user)
        patch_getSubset.assert_called_once_with(exp_id, 'named subset', user)

    @patch('ptmscout.views.dataset.dataset_explorer_view.compute_subset_enrichment')
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.annotations.getSubsetById')
    def test_fetch_subset_should_get_subset_by_id_if_id_given(self, patch_getSubset, patch_getExp, patch_compute):
        exp_id = 11
        exp = createMockExperiment(eid=11)
        user = createMockUser()
        subset = createMockSubset(name='named subset', owner_id=user.id, experiment_id=exp_id, foreground='some exp', background='some other exp')

        patch_getSubset.return_value = subset
        patch_getExp.return_value = exp
        subset.annotation_set_id = 10
        
        fetch_query = {
                       'experiment': str(exp_id),
                       'id': '100000',
                       'annotation_set_id': 10,
                       'motif_length': 5
                       }
                    
        request = DummyRequest()
        request.user = user
        
        request.body = json.dumps(fetch_query)
        
        patch_compute.return_value = {'some': 'enrichment analysis'}
        
        result = dataset_explorer_view.fetch_subset(request)
        
        self.assertEqual({'id':subset.id, 'some': 'enrichment analysis'}, result)
        patch_compute.assert_called_once_with(request, 10, exp, user, 'named subset', 'some exp', 'some other exp', 5)
        patch_getExp.assert_called_once_with(exp_id, user)
        patch_getSubset.assert_called_once_with(100000, exp_id)

        
        
class IntegrationTestDatasetExplorerView(IntegrationTestCase):
    def test_view_enrichment_page(self):
        result = self.ptmscoutapp.get("/experiments/28/subsets", status=200)
        result.mustcontain('Subsets')
        
    def perform_subset_fetch_test(self, exp_id, foreground_query, background_query = 'experiment'):
        query_expression = {
                            'experiment': exp_id,
                            'name': 'Subset 1',
                            'background': background_query,
                            'foreground': foreground_query
                           }
        
        result = self.ptmscoutapp.post_json("/webservice/subsets/query", params=query_expression, status=200)
        return result.json

    def test_subset_fetch_with_pfam(self):
        foreground_query = [
                            ['nop', ['metadata', 'Pfam-Site', 'eq', 'Pkinase']],
                            ['and', ['metadata', 'Pfam-Domain', 'eq', 'Pkinase']]
                            ]
        
        result = self.perform_subset_fetch_test(28, foreground_query)
        
        self.assertEqual(8, result['foreground']['peptides'])      
        self.assertEqual(7, result['foreground']['proteins'])      
        self.assertEqual(8, result['foreground']['sites'])
        
    def test_subset_fetch_integration(self):
        foreground_query = [
               ['nop', [ 'metadata', 'GO-Biological Process', 'eq', 'GO:0000187: activation of MAPK activity' ]]
            ]
        
        result = self.perform_subset_fetch_test(28, foreground_query)
                
        self.assertEqual(foreground_query, result['foreground']['query'])
        self.assertEqual(4, result['foreground']['proteins'])
        self.assertEqual(7, result['foreground']['peptides'])
        self.assertEqual(7, result['foreground']['sites'])
        
        self.assertEqual('experiment', result['background']['query'])
        self.assertEqual(52, result['background']['proteins'])
        self.assertEqual(68, result['background']['peptides'])
        self.assertEqual(70, result['background']['sites'])
        
        self.assertEqual('Subset 1', result['name'])
        self.assertEqual(28, result['experiment'])

        exp_peptides = [
                        [399228, u'NP_002737', u'MAPK3', u'http://localhost/proteins/35329?experiment_id=28', u'HTGFLTEyVATRWYR', u'Y204', u'Phosphotyrosine', {u'kinase': True, u'loop': u'K', 'mutant':False}], 
                        [399229, u'NP_002736', u'MAPK1', u'http://localhost/proteins/35548?experiment_id=28', u'HTGFLTEyVATRWYR', u'Y187', u'Phosphotyrosine', {u'kinase': True, u'loop': u'K', 'mutant':False}], 
                        [399230, u'NP_002736', u'MAPK1', u'http://localhost/proteins/35548?experiment_id=28', u'HDHTGFLtEYVATRW', u'T185', u'Phosphothreonine', {u'kinase': True, u'loop': u'K', 'mutant':False}], 
                        [399230, u'NP_002736', u'MAPK1', u'http://localhost/proteins/35548?experiment_id=28', u'HTGFLTEyVATRWYR', u'Y187', u'Phosphotyrosine', {u'kinase': True, u'loop': u'K', 'mutant':False}], 
                        [399247, u'Q16539', u'MAPK14', u'http://localhost/proteins/13914?experiment_id=28', u'TDDEMTGyVATRWYR', u'Y182', u'Phosphotyrosine', {u'kinase': True, u'loop': u'K', 'mutant':False}], 
                        [399263, u'NP_003020', u'SHC1', u'http://localhost/proteins/35338?experiment_id=28', u'EEPPDHQyYNDFPGK', u'Y239', u'Phosphotyrosine', {u'kinase': False, u'loop': None, 'mutant':False}], 
                        [399264, u'NP_003020', u'SHC1', u'http://localhost/proteins/35338?experiment_id=28', u'ELFDDPSyVNVQNLD', u'Y318', u'Phosphotyrosine', {u'kinase': False, u'loop': None, 'mutant':False}], 
                        [399265, u'NP_003020', u'SHC1', u'http://localhost/proteins/35338?experiment_id=28', u'EEPPDHQyYNDFPGK', u'Y239', u'Phosphotyrosine', {u'kinase': False, u'loop': None, 'mutant':False}], 
                        [399265, u'NP_003020', u'SHC1', u'http://localhost/proteins/35338?experiment_id=28', u'EPPDHQYyNDFPGKE', u'Y240', u'Phosphotyrosine', {u'kinase': False, u'loop': None, 'mutant':False}]
                    ]
        self.assertEqual(exp_peptides, result['peptides'])
        
    def test_subset_save_integration(self):
        from ptmscout.database import DBSession
        exp_id = 28
        foreground_query = [
               ['nop', [ 'metadata', 'GO-Biological Process', 'eq', 'GO:0000187: activation of MAPK activity' ]]
            ]
        background_query = 'experiment'
        
        query_expression = {
                            'experiment': exp_id,
                            'name': 'BP:MAPK',
                            'background': background_query,
                            'foreground': foreground_query}        
        
        result = self.ptmscoutapp.post_json("/webservice/subsets/save", params=query_expression, status=200)
        result = result.json
        
        
        DBSession.flush()
        
        subset = annotations.getSubsetByName(28, 'BP:MAPK', self.bot.user)
        
        self.assertEqual(query_expression['foreground'], subset.foreground_query)
        self.assertEqual(query_expression['background'], subset.background_query)
        
        foreground_query = [
                            ['nop', [ 'metadata', 'GO-Biological Process', 'eq', 'GO:0000186: activation of MAPKK activity' ]]
                            ]
        background_query = [
                            ['nop', ['subset', 'in', 'BP:MAPK']]
                            ]
        
        query_expression = {
                            'experiment':exp_id,
                            'name': 'BP:MAPK',
                            'background': background_query,
                            'foreground': foreground_query
                            }
        
        result = self.ptmscoutapp.post_json("/webservice/subsets/query", params = query_expression, status=200)
        result = result.json
        
        self.assertEqual(3, result['foreground']['peptides'])
        self.assertEqual(7, result['background']['peptides'])
