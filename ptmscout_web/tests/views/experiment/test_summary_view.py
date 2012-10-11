import unittest
from pyramid import testing
from ptmscout.views.experiment.summary_view import experiment_summary_view,\
    summarize_measurements, create_sequence_profile
from pyramid.testing import DummyRequest
from ptmscout.config import strings
from tests.views.mocking import createMockExperiment, createMockProtein,\
    createMockMeasurement, createMockPhosphopep
from mock import patch
import json
import base64
import math
from tests.PTMScoutTestCase import IntegrationTestCase

class SummaryViewIntegrationTests(IntegrationTestCase):
    def test_view_integration(self):
        self.ptmscoutapp.get('/experiments/26/summary')
        self.ptmscoutapp.get('/experiments/26')

class SummaryViewsTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        
    def tearDown(self):
        testing.tearDown()
    
    def test_get_sequence_profile(self):
        p1 = createMockProtein()
        
        m1 = createMockMeasurement(p1.id, 1)
        m2 = createMockMeasurement(p1.id, 1)
        
        pep1 = createMockPhosphopep(p1.id)
        pep2 = createMockPhosphopep(p1.id)
        pep3 = createMockPhosphopep(p1.id)
        m1.phosphopeps.append(pep1)
        m1.phosphopeps.append(pep2)
        m2.phosphopeps.append(pep3)
        
        pep1.pep_aligned='LKKVVALyDYMPMNA'
        pep2.pep_aligned='VSHWQQQsYLDSGIH'
        pep3.pep_aligned='VATWTAQsLLGSGIP'
        
        mods = [m1,m2]
        result = create_sequence_profile(mods)
        
        en = 19 / (2 * math.log(2) * 3)
        R21 = math.log(20, 2) + ( (2/3.0) * math.log((2/3.0), 2) + (1/3.0) * math.log((1/3.0), 2) ) - en
        R111 = math.log(20, 2) + ( (1/3.0) * math.log((1/3.0), 2) + (1/3.0) * math.log((1/3.0), 2) + (1/3.0) * math.log((1/3.0), 2) ) - en
        
        expected_peps = []
        expected_peps.append({'R':R21, 'f':[('V', 2, 0),('L', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('A', 1, 0),('K', 1, 1),('S', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('H', 1, 0),('K', 1, 1),('T', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('W', 2, 0),('V', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('Q', 1, 0),('T', 1, 1),('V', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('A', 2, 0),('Q', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('Q', 2, 0),('L', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('S', 2, 0),('Y', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('Y', 1, 0),('D', 1, 1),('L',1, 2)]})
        expected_peps.append({'R':R21, 'f':[('L', 2, 0),('Y', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('M', 1, 0),('D', 1, 1),('G',1, 2)]})
        expected_peps.append({'R':R21, 'f':[('S', 2, 0),('P', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('G', 2, 0),('M', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('I', 2, 0),('N', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('A', 1, 0),('H', 1, 1),('P',1, 2)]})
        
        expected = {'total': 3, 'frequencies': expected_peps}
        
        self.assertEqual(expected['total'], result['total'])
        self.assertEqual( len(expected['frequencies']), len(result['frequencies']) )
        
        for i in xrange(0, 15):
            exp = expected['frequencies'][i]
            res = result['frequencies'][i]
            
            self.assertAlmostEqual(exp['R'], res['R'], 6)
            self.assertEqual(exp['f'], res['f'])        
        
    
    def test_summarize_measurements(self):
        p1 = createMockProtein()
        m1 = createMockMeasurement(p1.id, 1)
        m2 = createMockMeasurement(p1.id, 1)
        
        m1.protein = p1
        m2.protein = p1
        
        pep1 = createMockPhosphopep(p1.id)
        pep2 = createMockPhosphopep(p1.id)
        pep3 = createMockPhosphopep(p1.id)
        
        pep1.site_type='Y'
        pep2.site_type='T'
        pep3.site_type='S'
        
        m1.phosphopeps.append(pep1)
        m1.phosphopeps.append(pep2)
        m2.phosphopeps.append(pep3)
        
        mods = [m1,m2]
        result = summarize_measurements(mods)
        
        self.assertEqual(3, result['modifications'])
        self.assertEqual(2, result['measured'])
        self.assertEqual(1, result['proteins'])
        
        self.assertEqual({'S':1, 'T':1, 'Y':1}, result['by_residue'])
        self.assertEqual({'homo sapiens': 3}, result['by_species'])
        self.assertEqual({'phosphorylation':3}, result['by_type'])

    @patch('ptmscout.views.experiment.summary_view.create_sequence_profile')
    @patch('ptmscout.views.experiment.summary_view.summarize_measurements')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_summary_view(self, patch_getExperiment, patch_getMeasurements, patch_measurementSummary, patch_sequenceProfile):
        request = DummyRequest()
        exp_id = 1
        request.matchdict['id'] = exp_id
        request.user = None
        
        sequence_profile = ["some list of items"]
        patch_sequenceProfile.return_value = sequence_profile
        
        measurements = [1,2,3,4]
        patch_getMeasurements.return_value = measurements
        
        mock_measurement_summary = {"a pep":"type dictionary"}
        patch_measurementSummary.return_value = mock_measurement_summary
        
        exp = createMockExperiment(exp_id, 0, 0)
        patch_getExperiment.return_value = exp
        
        result = experiment_summary_view(request)

        patch_getExperiment.assert_called_once_with(exp_id, None)
        patch_getMeasurements.assert_called_once_with(exp_id, request.user)
        
        patch_sequenceProfile.assert_called_once_with(measurements)        
        patch_measurementSummary.assert_called_once_with(measurements)
        
        self.assertEqual(exp, result['experiment'])

        self.assertEqual(mock_measurement_summary, result['measurement_summary'])
        self.assertEqual(base64.b64encode(json.dumps(sequence_profile)), result['sequence_profile'])
        
        self.assertEqual(0, result['error_count'])
        self.assertEqual(0, result['rejected_peptides'])
        self.assertEqual(strings.experiment_summary_page_title % (exp.name), result['pageTitle'])