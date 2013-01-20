from mock import patch
from ptmscout.config import strings
from ptmscout.database.taxonomies import Species
from ptmscout.views.protein.search_view import protein_search_view
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockProtein, createMockMeasurement, \
    createMockPeptide, createMockUser, createMockPTM,\
    createMockPeptideModification
from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase
import urllib

class TestProteinSearchViewIntegration(IntegrationTestCase):
#    def setUp(self):
#        IntegrationTestCase.setUp(self, sql_echo=True)

    def test_integration(self):
        result = self.ptmscoutapp.get('/proteins', {'acc_search':"ACK1", 'submitted':"true", 'species':"all"}, status=200)
        result.mustcontain("Showing 1 - 4 of 4 results")

    def test_integration_pep_search(self):
        result = self.ptmscoutapp.get('/proteins', {'pep_search':"RGR", 'submitted':"true", 'species':"all", 'page':'4'}, status=200)
        result.mustcontain("Showing 151 - 200 of 3585 results")

class TestProteinSearchViews(UnitTestCase):
       
    @patch('ptmscout.database.taxonomies.getAllSpecies')
    def test_protein_search_view_should_get_species_list_and_show_search_forms(self, patch_getSpecies):
        request = DummyRequest()
        
        species_list = [Species('mus musculus'), Species('homo sapiens')]
        patch_getSpecies.return_value = species_list
        
        result = protein_search_view(request)
        
        patch_getSpecies.assert_called_with()
        self.assertEqual(strings.protein_search_page_title, result['pageTitle'])
        self.assertEqual([], result['proteins'])
        self.assertEqual(False, result['submitted'])
        
    @patch('ptmscout.database.taxonomies.getAllSpecies')
    def test_protein_search_view_should_not_submit_search_if_protein_field_empty(self, patch_getSpecies):
        request = DummyRequest()
        request.GET['submitted'] = "true"
        request.GET['acc_search'] = ""
        request.GET['pep_search'] = ""
        request.GET['stringency'] = "2"
        request.GET['species'] = "homo sapiens"
        
        species_list = [Species('mus musculus'), Species('homo sapiens')]
        patch_getSpecies.return_value = species_list
        
        result = protein_search_view(request)
        
        patch_getSpecies.assert_called_with()
        self.assertEqual(strings.protein_search_page_title, result['pageTitle'])
        self.assertEqual([], result['proteins'])
        self.assertEqual(True, result['submitted'])
    

    def run_test_of_search(self, request, patch_getProteins, patch_getSpecies, patch_getMods):
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
        patch_getProteins.return_value = ( len(protein_list), protein_list )
        patch_getMods.return_value = mod_list
        
        result = protein_search_view(request)
        
        patch_getMods.assert_called_with(p1.id, ptm_user)
        patch_getProteins.assert_called_with(exp_id=None, search=request.GET['acc_search'],
            species=request.GET['species'], sequence=request.GET['pep_search']
        , page=(50,0))
        patch_getSpecies.assert_called_with()
        self.assertEqual(strings.protein_search_page_title, result['pageTitle'])
        
        sorted_peps = sorted([pep1, pep2, pep3], key=lambda pep: pep.getName())
        parsed_peps = [{'site':pep.getName(), 'peptide':pep.getPeptide()} for pep in sorted_peps]
        
        self.assertEqual(protein_list, result['proteins'])
        self.assertEqual(True, result['submitted'])

    @patch('ptmscout.database.modifications.getMeasuredPeptidesByProtein')
    @patch('ptmscout.database.protein.searchProteins')
    @patch('ptmscout.database.taxonomies.getAllSpecies')
    def test_protein_search_view_should_process_search_form(self, patch_getSpecies, patch_getProteins, patch_getMods):
        request = DummyRequest()
        request.GET['submitted'] = "true"
        request.GET['acc_search'] = "ACK1"
        request.GET['pep_search'] = ""
        request.GET['stringency'] = "2"
        request.GET['species'] = "homo sapiens"
        
        self.run_test_of_search(request, patch_getProteins, patch_getSpecies, patch_getMods)


    @patch('ptmscout.database.modifications.getMeasuredPeptidesByProtein')
    @patch('ptmscout.database.protein.searchProteins')
    @patch('ptmscout.database.taxonomies.getAllSpecies')
    def test_protein_search_view_should_process_search_form(self, patch_getSpecies, patch_getProteins, patch_getMods):
        request = DummyRequest()
        request.GET['submitted'] = "true"
        request.GET['acc_search'] = ""
        request.GET['pep_search'] = "ASDFKJEC"
        request.GET['stringency'] = "2"
        request.GET['species'] = "homo sapiens"
        
        self.run_test_of_search(request, patch_getProteins, patch_getSpecies, patch_getMods)
