from pyramid.testing import DummyRequest
from tests.views.mocking import createMockUser, createMockExperiment,\
    createMockError
from ptmscout.views.experiment.errors_view import experiment_errors_view
from mock import patch
from ptmscout.config import strings
from ptmscout.database import experiment
from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.httpexceptions import HTTPForbidden

class TestExperimentErrorsView(UnitTestCase):
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_view_should_raise_forbidden_if_user_not_owner(self, patch_getExp):
        user = createMockUser()
        
        exp = createMockExperiment()
        request = DummyRequest()
        request.matchdict['id'] = str(exp.id)
        request.user = user
        
        createMockError(line=4, message="Error 2", accession='acc3', peptide='pep4', experiment=exp)
        createMockError(line=1, message="Error 1", accession='acc', peptide='pep2', experiment=exp)
        
        user.experimentOwner.return_value = False
        patch_getExp.return_value = exp
        
        try:
            experiment_errors_view(request)
        except HTTPForbidden:
            pass
        else:
            self.fail("Expected exception: HTTPForbidden")
        
        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_view_should_get_error_listing_and_user_owner_status(self, patch_getExp):
        user = createMockUser()
        
        exp = createMockExperiment()
        request = DummyRequest()
        request.matchdict['id'] = str(exp.id)
        request.user = user
        
        createMockError(line=4, message="Error 2", accession='acc3', peptide='pep4', experiment=exp)
        createMockError(line=1, message="Error 1", accession='acc', peptide='pep2', experiment=exp)
        
        user.experimentOwner.return_value = True
        patch_getExp.return_value = exp
                        
        result = experiment_errors_view(request)
        
        patch_getExp.assert_called_once_with(exp.id, user)
        user.experimentOwner.assert_called_once_with(exp)
        
        self.assertEqual(True, result['user_owner'])
        
        expected_errors = [(1, 'acc', 'pep2',"Error 1"),
                           (4, 'acc3', 'pep4',"Error 2")]
        
        self.assertEqual(strings.experiment_errors_page_title % (exp.name), result['pageTitle'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(expected_errors,result['errors'])
        self.assertEqual(2,result['error_count'])
        self.assertEqual(2,result['rejected_peptides'])
        
        
class IntegrationTestExperimentErrorsView(IntegrationTestCase):
    def test_view_integration(self):
        self.bot.acquire_experiments([28])
        experiment.createExperimentError(28, 3, 'acc', 'pep', 'message')
        result = self.ptmscoutapp.get("/experiments/28/errors", status=200)
        
        result.mustcontain('http://localhost/upload?load_type=append&parent_experiment=28')
        
    def test_view_integration_with_datasets(self):
        self.bot.acquire_experiments([28])
        
        exp = experiment.getExperimentById(28, self.bot.user)
        exp.type='dataset'
        exp.saveExperiment()
        
        experiment.createExperimentError(28, 3, 'acc', 'pep', 'message')
        result = self.ptmscoutapp.get("/experiments/28/errors", status=200)
        
        result.mustcontain('http://localhost/dataset/upload?load_type=append&parent_experiment=28')
