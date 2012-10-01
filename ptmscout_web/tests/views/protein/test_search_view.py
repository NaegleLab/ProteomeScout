from mock import patch
from ptmscout.config import strings
from ptmscout.database.protein import Species
from ptmscout.views.protein.search_view import protein_search_view
from pyramid import testing
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockProtein, createMockModification, \
    createMockPhosphopep, createMockUser
import unittest

class TestProteinSearchViews(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
    @patch('ptmscout.database.protein.getAllSpecies')
    def test_protein_search_view_should_get_species_list_and_show_search_forms(self, patch_getSpecies):
        request = DummyRequest()
        
        species_list = [Species('mus musculus'), Species('homo sapiens')]
        patch_getSpecies.return_value = species_list
        
        result = protein_search_view(request)
        
        patch_getSpecies.assert_called_with()
        self.assertEqual(strings.protein_search_page_title, result['pageTitle'])
        self.assertEqual(['mus musculus', 'homo sapiens'], result['species_list'])
        self.assertEqual(None, result['selected_species'])
        self.assertEqual([], result['proteins'])
        self.assertEqual(False, result['include_predictions'])
        self.assertEqual({}, result['modifications'])
        self.assertEqual(False, result['submitted'])
        
    @patch('ptmscout.database.protein.getAllSpecies')
    def test_protein_search_view_should_not_submit_search_if_protein_field_empty(self, patch_getSpecies):
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.POST['acc_search'] = ""
        request.POST['stringency'] = "2"
        request.POST['species'] = "homo sapiens"
        
        species_list = [Species('mus musculus'), Species('homo sapiens')]
        patch_getSpecies.return_value = species_list
        
        result = protein_search_view(request)
        
        patch_getSpecies.assert_called_with()
        self.assertEqual(strings.protein_search_page_title, result['pageTitle'])
        self.assertEqual(['mus musculus', 'homo sapiens'], result['species_list'])
        
        self.assertEqual("", result['acc_search'])
        self.assertEqual("2", result['stringency'])
        self.assertEqual("homo sapiens", result['selected_species'])
        self.assertEqual([], result['proteins'])
        self.assertEqual(False, result['include_predictions'])
        self.assertEqual({}, result['modifications'])
        self.assertEqual(True, result['submitted'])
    
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByProtein')
    @patch('ptmscout.database.protein.getProteinsByAccession')
    @patch('ptmscout.database.protein.getAllSpecies')
    def test_protein_search_view_should_process_search_form(self, patch_getSpecies, patch_getProteins, patch_getMods):
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.POST['acc_search'] = "ACK1"
        request.POST['stringency'] = "2"
        request.POST['species'] = "homo sapiens"
        
        species_list = [Species('mus musculus'), Species('homo sapiens')]
        p1 = createMockProtein()
        protein_list = [p1]
        m1 = createMockModification(p1.id, 1)
        m2 = createMockModification(p1.id, 4)
        
        pep1 = createMockPhosphopep(p1.id)
        pep2 = createMockPhosphopep(p1.id)
        pep3 = createMockPhosphopep(p1.id)
        
        m1.phosphopeps.append(pep1)
        m1.phosphopeps.append(pep2)
        m2.phosphopeps.append(pep3)
        
        mod_list = [m1, m2]
        
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        patch_getSpecies.return_value = species_list
        patch_getProteins.return_value = protein_list
        patch_getMods.return_value = mod_list
        
        result = protein_search_view(request)
        
        patch_getMods.assert_called_with(p1.id, ptm_user)
        patch_getProteins.assert_called_with(["ACK1"], species="homo sapiens")
        patch_getSpecies.assert_called_with()
        self.assertEqual(strings.protein_search_page_title, result['pageTitle'])
        self.assertEqual(['mus musculus', 'homo sapiens'], result['species_list'])
        
        self.assertEqual("ACK1", result['acc_search'])
        self.assertEqual("2", result['stringency'])
        self.assertEqual("homo sapiens", result['selected_species'])
        
        sorted_peps = sorted([pep1, pep2, pep3], key=lambda pep: pep.getName())
        parsed_peps = [{'site':pep.getName(), 'peptide':pep.getPeptide()} for pep in sorted_peps]
        
        self.assertEqual(protein_list, result['proteins'])
        self.assertEqual(False, result['include_predictions'])
        self.assertEqual({p1.id: parsed_peps}, result['modifications'])
        self.assertEqual(True, result['submitted'])