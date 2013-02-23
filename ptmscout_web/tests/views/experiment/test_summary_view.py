from ptmscout.views.experiment.summary_view import experiment_summary_view,\
    summarize_measurements
from pyramid.testing import DummyRequest
from ptmscout.config import strings
from tests.views.mocking import createMockExperiment, createMockProtein,\
    createMockMeasurement, createMockPeptide, createMockError,\
    createMockPeptideModification, createMockPTM, createMockUser
from mock import patch
import json
import base64
from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase

class SummaryViewIntegrationTests(IntegrationTestCase):
    def test_view_integration(self):
        self.ptmscoutapp.get('/experiments/26/summary')
        self.ptmscoutapp.get('/experiments/26')

class SummaryViewsTests(UnitTestCase):
    def test_summarize_measurements(self):
        p1 = createMockProtein()
        m1 = createMockMeasurement(p1.id, 1)
        m2 = createMockMeasurement(p1.id, 1)
        
        m1.protein = p1
        m2.protein = p1
        
        pep1 = createMockPeptide(p1.id)
        pep2 = createMockPeptide(p1.id)
        pep3 = createMockPeptide(p1.id)
        
        pep1.site_type='Y'
        pep2.site_type='T'
        pep3.site_type='S'
        
        mod = createMockPTM(name="phosphorylation")
        mod2 = createMockPTM(name="phosphoserine",parent=mod)
        
        createMockPeptideModification(m1, pep1, mod)
        createMockPeptideModification(m1, pep2, mod)
        createMockPeptideModification(m2, pep3, mod2)
        
        mods = [m1,m2]
        result = summarize_measurements(mods)
        
        self.assertEqual(3, result['modifications'])
        self.assertEqual(2, result['measured'])
        self.assertEqual(1, result['proteins'])
        
        self.assertEqual({'S':1, 'T':1, 'Y':1}, result['by_residue'])
        self.assertEqual({'homo sapiens': 3}, result['by_species'])
        self.assertEqual({'phosphorylation':3}, result['by_type'])

    @patch('ptmscout.utils.protein_utils.create_sequence_profile')
    @patch('ptmscout.views.experiment.summary_view.summarize_measurements')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_summary_view(self, patch_getExperiment, patch_getMeasurements, patch_measurementSummary, patch_sequenceProfile):
        request = DummyRequest()
        exp_id = 1
        request.matchdict['id'] = exp_id
        request.user = createMockUser()
        
        sequence_profile = ["some list of items"]
        patch_sequenceProfile.return_value = sequence_profile
        
        measurements = [1,2,3,4]
        patch_getMeasurements.return_value = measurements
        
        mock_measurement_summary = {"a pep":"type dictionary"}
        patch_measurementSummary.return_value = mock_measurement_summary
        
        exp = createMockExperiment(exp_id, 0, 0)
        
        createMockError(1, "Not working", accession="QUP123", peptide="adsfgdfgt", experiment=exp)
        createMockError(2, "Not working", accession="QUP123", peptide="adsfgdfgt", experiment=exp)
        createMockError(3, "Not working", accession="KPU321", peptide="wertsdfgf", experiment=exp)
        
        patch_getExperiment.return_value = exp
        
        request.user.experimentOwner.return_value = True
        
        result = experiment_summary_view(request)

        request.user.experimentOwner.assert_called_once_with(exp)

        patch_getExperiment.assert_called_once_with(exp_id, request.user)
        patch_getMeasurements.assert_called_once_with(exp_id, request.user)
        
        patch_sequenceProfile.assert_called_once_with(measurements)        
        patch_measurementSummary.assert_called_once_with(measurements)
        
        self.assertEqual(exp, result['experiment'])

        self.assertEqual(mock_measurement_summary, result['measurement_summary'])
        self.assertEqual(base64.b64encode(json.dumps(sequence_profile)), result['sequence_profile'])
        
        self.assertEqual(True, result['user_owner'])
        self.assertEqual(2, result['rejected_peptides'])
        self.assertEqual(strings.experiment_summary_page_title % (exp.name), result['pageTitle'])
