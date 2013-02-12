from mock import patch, call
from ptmscout.config import strings
from ptmscout.database import Species
from ptmscout.views.experiment.browse_view import browse_experiment
from pyramid import testing
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockExperiment, createMockProtein, \
    createMockMeasurement, createMockPeptide, createMockScansite,\
    createMockUser, createMockPTM, createMockPeptideModification
import unittest
from tests.PTMScoutTestCase import IntegrationTestCase

class ExperimentBrowseViewIntegrationTests(IntegrationTestCase):
    def test_experiment_browse_view(self):
        result = self.ptmscoutapp.get('/experiments/26/browse')
        result.mustcontain('Showing 1 - 50 of 58')
        
    def test_experiment_browse_search(self):
        result = self.ptmscoutapp.get('/experiments/26/browse', {'submitted':"true", 'acc_search':"ACK1"})
        result.mustcontain('Showing 1 - 1 of 1 results')

    def test_experiment_browse_search_none(self):
        result = self.ptmscoutapp.get('/experiments/26/browse', {'submitted':"true", 'acc_search':"", 'pep_search':""})
        result.mustcontain('Showing 0 - 0 of 0 results')
        result.mustcontain('At least one of Protein, Peptide are required')



class ExperimentBrowseViewTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
    def setup_experiment_browse_test(self, patch_getExperiment, patch_getProteins, patch_getMods):
        mock_experiment = createMockExperiment(1, 0, 0)
        patch_getExperiment.return_value = mock_experiment

        mock_experiment.name = "experiment name"
        
        request = DummyRequest()
        
        ptm_user = createMockUser()
        request.user = ptm_user
        
        p1 = createMockProtein()
        p2 = createMockProtein()
        
        protein_list = sorted([p1, p2], key=lambda item: item.acc_gene)
        p1,p2 = protein_list
        patch_getProteins.return_value = len(protein_list), protein_list

        m1 = createMockMeasurement(p1.id, 1)
        m2 = createMockMeasurement(p1.id, 4)
        m3 = createMockMeasurement(p2.id, 6)
        m1.protein = p1
        m2.protein = p1
        m3.protein = p2
        m1.experiment_id = mock_experiment.id
        m2.experiment_id = mock_experiment.id
        m3.experiment_id = mock_experiment.id
        
        pep1 = createMockPeptide(p1.id)
        pep2 = createMockPeptide(p1.id)
        pep3 = createMockPeptide(p1.id)
        pep4 = createMockPeptide(p2.id)
        
        mod = createMockPTM()
        
        createMockPeptideModification(m1, pep1, mod)
        createMockPeptideModification(m1, pep2, mod)
        createMockPeptideModification(m2, pep3, mod)
        createMockPeptideModification(m3, pep4, mod)
        
        s1 = createMockScansite(pep1.id)
        s2 = createMockScansite(pep2.id)
        s3 = createMockScansite(pep3.id)
        
        s1.source='scansite'
        s2.source='scansite_kinase'
        s3.source='scansite'
        
        pep1.predictions = [s1, s2]
        pep2.predictions = [s1]
        pep4.predictions = [s2, s3]
        
        mod_list = [m1, m2, m3]
        patch_getMods.return_value = mod_list
        
        return request, p1, p2, mock_experiment, pep1, pep2, pep3, pep4, s1, s2, s3, protein_list

    @patch('ptmscout.database.taxonomies.getAllSpecies')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByProtein')
    @patch('ptmscout.database.protein.searchProteins')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_browse_view_when_searching(self, patch_getExperiment, patch_getProteins, patch_getMods, patch_getSpecies):
        species_list = [Species('mus musculus'), Species('homo sapiens')]
        patch_getSpecies.return_value = species_list
        request, p1, p2, mock_experiment, pep1, pep2, pep3, pep4, _, s2, s3, protein_list = self.setup_experiment_browse_test(patch_getExperiment, patch_getProteins, patch_getMods)
        
        request.GET['submitted'] = "true"
        request.GET['acc_search'] = "ACK1"
        request.matchdict['id'] = 1
        
        result = browse_experiment(request)
        
        patch_getProteins.assert_called_once_with(exp_id=1, search='ACK1', page=(50, 0), species=None, sequence='')
        patch_getExperiment.assert_called_once_with(1, request.user)
        self.assertIn(call(p1.id, request.user), patch_getMods.call_args_list)
        self.assertIn(call(p2.id, request.user), patch_getMods.call_args_list)

        self.assertEqual(mock_experiment, result['experiment'])
        self.assertEqual(strings.experiment_browse_page_title % ("experiment name"), result['pageTitle'])
        
        self.assertEqual(protein_list, result['proteins'])
        self.assertEqual(True, result['submitted'])
        
    @patch('ptmscout.database.taxonomies.getAllSpecies')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByProtein')
    @patch('ptmscout.database.protein.getProteinsByExperiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_browse_view_when_no_search_parameters(self, patch_getExperiment, patch_getProteins, patch_getMods, patch_getSpecies):
        species_list = [Species('mus musculus'), Species('homo sapiens')]
        patch_getSpecies.return_value = species_list
        request, p1, p2, mock_experiment, pep1, pep2, pep3, pep4, _, s2, s3, protein_list = self.setup_experiment_browse_test(patch_getExperiment, patch_getProteins, patch_getMods)
        request.matchdict['id'] = 1
        
        result = browse_experiment(request)
        
        patch_getExperiment.assert_called_once_with(1, request.user)
        self.assertIn(call(p1.id, request.user), patch_getMods.call_args_list)
        self.assertIn(call(p2.id, request.user), patch_getMods.call_args_list)

        self.assertEqual(mock_experiment, result['experiment'])
        self.assertEqual(strings.experiment_browse_page_title % ("experiment name"), result['pageTitle'])
        
        self.assertEqual(protein_list, result['proteins'])
        self.assertEqual(False, result['submitted'])
