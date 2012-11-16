from mock import patch
from ptmscout.config import strings
from ptmscout.views.protein.GO_view import protein_gene_ontology_view
from pyramid import testing
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockProtein, createMockGO, createMockGOE
import unittest

from tests.PTMScoutTestCase import IntegrationTestCase

class ProteinGOIntegrationTests(IntegrationTestCase):
    def test_protein_GO(self):
        self.ptmscoutapp.get('/proteins/35546/GO')
        



class TestProteinGOViews(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
    @patch('ptmscout.database.protein.getProteinById')
    def test_protein_GO_view_should_get_gene_ontology_info_for_protein(self, patch_getProtein):
        request = DummyRequest()

        mock_prot = createMockProtein()
        request.matchdict['id'] = str(mock_prot.id)
        patch_getProtein.return_value = mock_prot
        
        
        mock_prot.GO_terms.append(createMockGOE(mock_prot.id, createMockGO('P')))
        mock_prot.GO_terms.append(createMockGOE(mock_prot.id, createMockGO('C')))
        mock_prot.GO_terms.append(createMockGOE(mock_prot.id, createMockGO('P')))
        mock_prot.GO_terms.append(createMockGOE(mock_prot.id, createMockGO('F')))
        mock_prot.GO_terms.append(createMockGOE(mock_prot.id, createMockGO('C')))
        mock_prot.GO_terms.append(createMockGOE(mock_prot.id, createMockGO('P')))
        mock_prot.GO_terms.append(createMockGOE(mock_prot.id, createMockGO('C')))
        mock_prot.GO_terms.append(createMockGOE(mock_prot.id, createMockGO('F')))
        
        result = protein_gene_ontology_view(request)
        
        patch_getProtein.assert_called_with(mock_prot.id)
        
        self.assertEqual(strings.protein_ontology_page_title, result['pageTitle'])
        self.assertEqual(mock_prot, result['protein'])
        self.assertEqual(2, len(result['F_terms']))
        self.assertEqual(3, len(result['P_terms']))
        self.assertEqual(3, len(result['C_terms']))