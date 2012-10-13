from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from tests.behave.steps.dataset_load_steps import set_form_defaults
from ptmscout.config import settings, strings
import os
from pyramid.testing import DummyRequest
from ptmscout.views.experiment.upload_status import upload_confirm_view
from tests.views.mocking import createMockExperiment, createMockUser
from mock import patch
from pyramid.httpexceptions import HTTPFound

class TestUploadStatusView(UnitTestCase):
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_start_upload_view_should_redirect_if_already_uploaded(self, patch_getExperiment):
        request = DummyRequest()
        request.matchdict['id'] = "26"
        request.user = createMockUser("username", "email", "password", 0) 
        
        exp = createMockExperiment(26, 0, 0)
        exp.ready = 1
        
        patch_getExperiment.return_value = exp
        
        f = upload_confirm_view(request)
        self.assertEqual(request.application_url + "/experiment/26", f.location)
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_start_upload_view_should_get_confirmation(self, patch_getExperiment):
        request = DummyRequest()
        request.matchdict['id'] = "26"
        request.user = createMockUser("username", "email", "password", 0) 
        
        exp = createMockExperiment(26, 0, 0)
        exp.ready = 0
        
        patch_getExperiment.return_value = exp
        
        result = upload_confirm_view(request)
        
        patch_getExperiment.assert_called_once_with(26, request.user, False)
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.experiment_upload_confirm_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_confirm_message, result['message'])
        self.assertEqual("false", result['confirm'])
        
        
        
class IntegrationTestUploadStatusView(IntegrationTestCase):
    def test_view_integration(self):
        self.bot.login()
        result = self.ptmscoutapp.get("/upload", status=200)
        
        set_form_defaults(result.form)
        
        result.form.set('load_type', "new")
        
        os.chdir(settings.ptmscout_path)
        filename = "tests/behave/data/datasetLoad_correctDataset.txt"
        f = open(filename, 'rb')
        filecontents = f.read()
        
        result.form.set('data_file', (filename, filecontents))
        result.form.set('description', "This is a correct dataset")
        
        result = result.form.submit().follow()
        
        result.mustcontain(strings.experiment_upload_confirm_message)
        
        result = result.form.submit()
        
        result.mustcontain(strings.experiment_upload_started_message % ("http://localhost/account/experiments"))