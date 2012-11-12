from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase
from mock import patch
from ptmscout.database import experiment
import os
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockExperiment, createMockUser,\
    createMockError
from ptmscout.views.experiment.download_view import download_experiment
from pyramid.httpexceptions import HTTPForbidden
from ptmscout.utils import uploadutils
from ptmscout.config import strings

class TestExperimentDownloadView(UnitTestCase):
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_download_experiment_should_raise_forbidden_to_non_owner(self, patch_getExperiment):
        data_file = 'test/test_dataset_formatted.txt'
        exp = createMockExperiment()
        exp.dataset = data_file
        
        createMockError(1, "1 Error", experiment=exp)
        createMockError(8, "2 Error", experiment=exp)
        createMockError(8, "3 Error", experiment=exp)
        createMockError(17, "4 Error", experiment=exp)
        
        user = createMockUser()
        user.myExperiments.return_value = []
        
        request = DummyRequest()
        request.matchdict['id'] = exp.id
        request.user = user
        
        patch_getExperiment.return_value = exp
        
        try:
            download_experiment(request)
        except HTTPForbidden:
            pass
        else:
            self.fail("Expected HTTPForbidden")
        
        patch_getExperiment.assert_called_once_with(exp.id, user)
        user.myExperiments.assert_called_once_with()
    
        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_download_experiment_view_result_should_contain_recorded_errors(self, patch_getExperiment):
        data_file = 'test/test_dataset_formatted.txt'
        exp = createMockExperiment()
        exp.dataset = data_file
        
        exp_headers, exp_data = uploadutils.load_header_and_data_rows(data_file, 100)
        exp_headers.insert(0, strings.experiment_upload_error_reasons_column_title)
        
        for i in xrange(0, len(exp_data)):
            exp_data[i].insert(0, "")
            
        exp_data[0][0] = "1 Error"
        exp_data[7][0] = "2 Error, 3 Error"
        exp_data[16][0] = "4 Error"
        
        exp_data = [row for row in exp_data if row[0] != ""]
        
        createMockError(1, "1 Error", experiment=exp)
        createMockError(8, "2 Error", experiment=exp)
        createMockError(8, "3 Error", experiment=exp)
        createMockError(17, "4 Error", experiment=exp)
        
        user = createMockUser()
        user.myExperiments.return_value = [exp]
        
        request = DummyRequest()
        request.GET['errors'] = "True"
        request.matchdict['id'] = exp.id
        request.user = user
        
        patch_getExperiment.return_value = exp
        
        result = download_experiment(request)
        
        patch_getExperiment.assert_called_once_with(exp.id, user)
        user.myExperiments.assert_called_once_with()
        
        
        self.assertEqual(exp_headers, result['header'])
        self.assertEqual(exp_data, result['data'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_download_experiment_view(self, patch_getExperiment):
        data_file = 'test/test_dataset_formatted.txt'
        exp = createMockExperiment()
        exp.dataset = data_file
        
        exp_headers, exp_data = uploadutils.load_header_and_data_rows(data_file, 100)
        
        createMockError(1, "1 Error", experiment=exp)
        createMockError(8, "2 Error", experiment=exp)
        createMockError(8, "3 Error", experiment=exp)
        createMockError(17, "4 Error", experiment=exp)
        
        user = createMockUser()
        user.myExperiments.return_value = [exp]
        
        request = DummyRequest()
        request.matchdict['id'] = exp.id
        request.user = user
        
        patch_getExperiment.return_value = exp
        
        result = download_experiment(request)
        
        patch_getExperiment.assert_called_once_with(exp.id, user)
        user.myExperiments.assert_called_once_with()
        
        self.assertEqual(exp_headers, result['header'])
        self.assertEqual(exp_data, result['data'])
        
        
class IntegrationTestExperimentDownloadView(IntegrationTestCase):
    def test_view_should_forbidden_if_not_logged_in(self):
        self.bot.logout()
        result = self.ptmscoutapp.get("/experiments/28/download", status=200)
        result.mustcontain("forbidden")
        
    def test_view_should_forbidden_if_not_owner(self):
        self.bot.login()
        result = self.ptmscoutapp.get("/experiments/28/download", status=200)
        result.mustcontain("forbidden")
        
    def test_view_should_succeed_if_owner(self):
        self.bot.login()
        
        exp = experiment.getExperimentById(28, secure=False)
        exp.dataset = os.path.join('test','test_dataset_formatted.txt')
        exp.saveExperiment()
        
        self.bot.acquire_experiments([28])
        
        self.ptmscoutapp.get("/experiments/28/download", status=200)
        