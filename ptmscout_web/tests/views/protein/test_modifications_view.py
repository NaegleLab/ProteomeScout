from mock import patch
from ptmscout.config import strings
from ptmscout.views.protein.modifications_view import protein_modifications_view
from pyramid import testing
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockProtein, createMockUser, \
    createMockMeasurement, createMockExperiment, createMockPeptide,\
    createMockPTM, createMockPeptideModification
import unittest
from tests.PTMScoutTestCase import IntegrationTestCase


class TestProteinModificationViewsIntegration(IntegrationTestCase):
    def test_integration(self):
        self.ptmscoutapp.get('/proteins/35546/modifications')

class TestProteinModificationViews(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByProtein')
    @patch('ptmscout.database.protein.getProteinById')
    def test_protein_modifications_view_should_show_modifications_for_protein(self, patch_getProtein, patch_getMods):
        request = DummyRequest()

        mock_prot = createMockProtein()
        request.matchdict['id'] = str(mock_prot.id)
        patch_getProtein.return_value = mock_prot
        
        mock_user = createMockUser("username", "email", "password", 1)
        request.user = mock_user
        
        mock_mod = createMockMeasurement(mock_prot.id, 24)
        mock_mod.experiment = createMockExperiment(2, 0, 0)
        
        mock_mod2 = createMockMeasurement(mock_prot.id, 25)
        mock_mod2.experiment = createMockExperiment(3, 0, 0)
        
        p1 = createMockPeptide(mock_prot.id)
        
        mod = createMockPTM()
        
        createMockPeptideModification(mock_mod, p1, mod)
        createMockPeptideModification(mock_mod2, p1, mod)
        
        p1.getName.return_value = "mock_pep_name"
        p1.getPeptide.return_value = "mock_peptide"
        
        mod_list = [mock_mod, mock_mod2]
        patch_getMods.return_value = mod_list 
        
        result = protein_modifications_view(request)
        
        patch_getProtein.assert_called_with(mock_prot.id)
        patch_getMods.assert_called_with(mock_prot.id, mock_user)
        
        self.assertEqual(strings.protein_modification_sites_page_title, result['pageTitle'])
        self.assertEqual(mock_prot, result['protein'])
        self.assertEqual([{'name':"mock_pep_name",'peptide':"mock_peptide",'experiments':set([mock_mod.experiment, mock_mod2.experiment])}], result['modification_sites'])
    