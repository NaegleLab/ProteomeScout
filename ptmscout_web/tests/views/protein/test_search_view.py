from mock import patch
from ptmscout.config import strings
from ptmscout.database.taxonomies import Species
from ptmscout.views.protein.search_view import protein_search_view
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockProtein, createMockMeasurement, \
    createMockPeptide, createMockUser, createMockPTM,\
    createMockPeptideModification
from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase

class TestProteinSearchViewIntegration(IntegrationTestCase):
    def test_integration(self):
                self.ptmscoutapp.get('/proteins', {'acc_search':"ACK1", 'stringency':"2", 'submitted':"1", 'species':"all"}, status=200)
                
class TestProteinSearchViews(UnitTestCase):
        
    @patch('ptmscout.database.taxonomies.getAllSpecies')
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
        
    @patch('ptmscout.database.taxonomies.getAllSpecies')
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
    @patch('ptmscout.database.protein.searchProteins')
    @patch('ptmscout.database.taxonomies.getAllSpecies')
    def test_protein_search_view_should_process_search_form(self, patch_getSpecies, patch_getProteins, patch_getMods):
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.POST['acc_search'] = "ACK1"
        request.POST['stringency'] = "2"
        request.POST['species'] = "homo sapiens"
        
        species_list = [Species('mus musculus'), Species('homo sapiens')]
        p1 = createMockProtein()
        protein_list = [p1]
        m1 = createMockMeasurement(p1.id, 1)
        m2 = createMockMeasurement(p1.id, 4)
        
        pep1 = createMockPeptide(p1.id)
        pep2 = createMockPeptide(p1.id)
        pep3 = createMockPeptide(p1.id)
        
        mod = createMockPTM()
        
        createMockPeptideModification( m1, pep1, mod )
        createMockPeptideModification( m1, pep2, mod )
        createMockPeptideModification( m2, pep3, mod )
        
        mod_list = [m1, m2]
        
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        patch_getSpecies.return_value = species_list
        patch_getProteins.return_value = len(protein_list), protein_list
        patch_getMods.return_value = mod_list
        
        result = protein_search_view(request)
        
        patch_getMods.assert_called_with(p1.id, ptm_user)
        patch_getProteins.assert_called_with("ACK1", species="homo sapiens")
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
