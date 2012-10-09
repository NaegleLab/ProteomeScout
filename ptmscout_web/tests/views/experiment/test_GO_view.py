from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
from mock import patch
from tests.views.mocking import createMockExperiment, createMockMeasurement,\
    createMockProtein, createMockGO
from ptmscout.views.experiment.GO_view import show_experiment_go_terms,\
    format_go_terms, build_go_annotation_tree
from ptmscout.config import strings

class TestExperimentGOView(UnitTestCase):
    
    def test_build_go_annotation_tree_should_build_json_tree(self):
        result = build_go_annotation_tree([])
        self.assertEqual("", result)
    
    def test_format_go_terms_should_build_tables(self):
        
        p1 = createMockProtein()
        p2 = createMockProtein()
        
        g1 = createMockGO('F')
        g2 = createMockGO('P')
        g3 = createMockGO('C')
        g4 = createMockGO('C')
        
        p1.GO_terms = [g1,g2,g3]
        p2.GO_terms = [g2,g4]
        
        m1 = createMockMeasurement(p1.id, 1)
        m2 = createMockMeasurement(p1.id, 1)
        m3 = createMockMeasurement(p2.id, 1)
        
        m1.protein = p1
        m2.protein = p2
        m3.protein = p2
        
        measures = [m1,m2,m3]
        
        result = format_go_terms(measures)
        
        self.assertEqual({'molecular_function':[(g1.GO, g1.term, 1)],
                          'biological_process':[(g2.GO, g2.term, 2)],
                          'cellular_component':[(g4.GO, g4.term, 1), (g3.GO, g3.term, 1)]}, result)
    
    @patch('ptmscout.views.experiment.GO_view.build_go_annotation_tree')
    @patch('ptmscout.views.experiment.GO_view.format_go_terms')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_view(self, patch_getExperiment, patch_getMeasurements, patch_format, patch_tree):
        request = DummyRequest()
        exp_id = 1
        request.matchdict['id'] = str(exp_id)
        request.user = None
        
        exp = createMockExperiment(exp_id,0,0)
        patch_getExperiment.return_value = exp
        
        measures = ["some","measurements"]
        patch_getMeasurements.return_value = measures
        
        formatted = ["some","formatted","measurements"]
        patch_format.return_value = formatted
        
        tree = "some_json_format_tree"
        patch_tree.return_value = tree
        
        result = show_experiment_go_terms(request)
        
        patch_format.assert_called_once_with(measures)
        patch_getMeasurements.assert_called_once_with(exp_id, request.user)
        patch_getExperiment.assert_called_once_with(exp_id, request.user)
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.experiment_GO_page_title % (exp.name), result['pageTitle'])
        self.assertEqual(formatted, result['go_tables'])
        self.assertEqual(tree, result['go_tree'])
        
        
        
        
class IntegrationTestExperimentGOView(IntegrationTestCase):
    def test_view_integration(self):
        self.ptmscoutapp.get("/experiments/26/GO", status=200)