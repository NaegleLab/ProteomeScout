import unittest
from pyramid.testing import DummyRequest
from ptmscout.proteins import protein_modifications_view, protein_search_view,\
    protein_expression_view, protein_gene_ontology_view,\
    protein_experiment_data_view
from ptmscout import strings
from tests.views.mocking import createMockProtein, createUserForTest,\
    createMockModification, createMockExperiment, createMockPhosphopep,\
    createMockProbe, createMockExpSample, createMockGO
from mock import patch
from ptmscout.database.protein import Species



class TestProteinViews(unittest.TestCase):
    
    @patch('ptmscout.database.protein.getProteinById')
    def test_protein_data_view_get_experiment_data_for_protein(self, patch_getProtein):
        request = DummyRequest()

        mock_prot = createMockProtein()
        request.matchdict['id'] = str(mock_prot.id)
        patch_getProtein.return_value = mock_prot
        
        result = protein_experiment_data_view(request)
        
        patch_getProtein.assert_called_with(mock_prot.id)
        
        self.assertEqual(strings.protein_data_page_title, result['pageTitle'])
        self.assertEqual(mock_prot, result['protein'])
    
    @patch('ptmscout.database.protein.getProteinById')
    def test_protein_GO_view_should_get_gene_ontology_info_for_protein(self, patch_getProtein):
        request = DummyRequest()

        mock_prot = createMockProtein()
        request.matchdict['id'] = str(mock_prot.id)
        patch_getProtein.return_value = mock_prot
        
        mock_prot.GO_terms.append(createMockGO('P'))
        mock_prot.GO_terms.append(createMockGO('C'))
        mock_prot.GO_terms.append(createMockGO('P'))
        mock_prot.GO_terms.append(createMockGO('F'))
        mock_prot.GO_terms.append(createMockGO('C'))
        mock_prot.GO_terms.append(createMockGO('P'))
        mock_prot.GO_terms.append(createMockGO('C'))
        mock_prot.GO_terms.append(createMockGO('F'))
        
        result = protein_gene_ontology_view(request)
        
        patch_getProtein.assert_called_with(mock_prot.id)
        
        self.assertEqual(strings.protein_ontology_page_title, result['pageTitle'])
        self.assertEqual(mock_prot, result['protein'])
        self.assertEqual(2, len(result['F_terms']))
        self.assertEqual(3, len(result['P_terms']))
        self.assertEqual(3, len(result['C_terms']))
        
    
    @patch('ptmscout.database.protein.getProteinById')
    def test_protein_expression_view_should_get_expression_data_for_protein_from_request_probeset(self, patch_getProtein):
        request = DummyRequest()

        mock_prot = createMockProtein()
        request.matchdict['id'] = str(mock_prot.id)
        patch_getProtein.return_value = mock_prot
        probe = createMockProbe()
        mock_prot.expression_probes = [probe]
        
        s1 = createMockExpSample(probe.id, 1, 'Human', 17, 'Tongue')
        s2 = createMockExpSample(probe.id, 1, 'Human', 1, '721 B Lymphocytes')
        s3 = createMockExpSample(probe.id, 3, 'NCI60', 20, 'Renal Carcinoma')
        s4 = createMockExpSample(probe.id, 3, 'NCI60', 23, 'Mammary')
        
        probe.samples = [s1,s2,s3,s4]
        
        result = protein_expression_view(request)
        
        patch_getProtein.assert_called_with(mock_prot.id)
        
        self.assertEqual(strings.protein_expression_page_title, result['pageTitle'])
        self.assertEqual(mock_prot, result['protein'])
        
        self.assertEqual([probe.probeset_id], result['probe_ids'])
        self.assertEqual(['Human', 'NCI60'], result['collections'])
        
        exp_data = [{'probeset':probe.probeset_id, 'collection':s.collection.name, 'tissue':s.tissue.name, 'value':s.value} for s in probe.samples]
        self.assertEqual(exp_data, result['expression_data'])
        
    
    @patch('ptmscout.database.modifications.getModificationsByProtein')
    @patch('ptmscout.database.protein.getProteinById')
    def test_protein_modifications_view_should_show_modifications_for_protein(self, patch_getProtein, patch_getMods):
        request = DummyRequest()

        mock_prot = createMockProtein()
        request.matchdict['id'] = str(mock_prot.id)
        patch_getProtein.return_value = mock_prot
        
        mock_user = createUserForTest("username", "email", "password", 1)
        request.user = mock_user
        
        mock_mod = createMockModification(mock_prot.id, 24)
        mock_mod.experiment = createMockExperiment(2, 0, 0)
        
        mock_mod2 = createMockModification(mock_prot.id, 25)
        mock_mod2.experiment = createMockExperiment(3, 0, 0)
        
        p1 = createMockPhosphopep(mock_prot.id)
        mock_mod.phosphopeps.append(p1)
        mock_mod2.phosphopeps.append(p1)
        
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
        self.assertEqual({}, result['modifications'])
        self.assertEqual(True, result['submitted'])
    
    @patch('ptmscout.database.modifications.getModificationsByProtein')
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
        
        ptm_user = createUserForTest("username", "email", "password", 1)
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
        self.assertEqual({p1.id: parsed_peps}, result['modifications'])
        self.assertEqual(True, result['submitted'])