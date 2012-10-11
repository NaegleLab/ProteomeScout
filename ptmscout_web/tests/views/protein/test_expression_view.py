from mock import patch
from ptmscout.config import strings
from ptmscout.views.protein.expression_view import protein_expression_view
from pyramid import testing
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockProtein, createMockProbe, \
    createMockExpSample
import unittest
from tests.PTMScoutTestCase import IntegrationTestCase

class TestProteinExpressionViewsIntegration(IntegrationTestCase):
    def test_integration(self):
        self.ptmscoutapp.get('/proteins/35546/expression')


class TestProteinExpressionViews(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
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
        