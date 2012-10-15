from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from tests.behave.steps.dataset_load_steps import set_form_defaults
from ptmscout.config import settings, strings
import os
from pyramid.testing import DummyRequest
from ptmscout.views.experiment.upload_status import upload_confirm_view,\
    UploadAlreadyStarted
from tests.views.mocking import createMockExperiment, createMockUser
from mock import patch
from pyramid.httpexceptions import HTTPFound

class TestUploadStatusView(UnitTestCase):
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_start_upload_view_should_redirect_if_already_started(self, patch_getExperiment):
        request = DummyRequest()
        request.matchdict['id'] = "26"
        request.user = createMockUser("username", "email", "password", 0) 
        
        exp = createMockExperiment(26, 0, 0, 'loading')
        
        patch_getExperiment.return_value = exp
        
        try:
            upload_confirm_view(request)
        except UploadAlreadyStarted, uas:
            self.assertEqual(26, uas.eid)
        else:
            self.fail("Expected exception UploadAlreadyStarted")
    
        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_start_upload_view_should_redirect_if_already_uploaded(self, patch_getExperiment):
        request = DummyRequest()
        request.matchdict['id'] = "26"
        request.user = createMockUser("username", "email", "password", 0) 
        
        exp = createMockExperiment(26, 0, 0, 'loaded')
        
        patch_getExperiment.return_value = exp
        
        f = upload_confirm_view(request)
        self.assertEqual(request.application_url + "/experiment/26", f.location)

    @patch('ptmscout.utils.data.start_upload')        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_start_upload_view_should_start_job_and_display_confirmation(self, patch_getExperiment, patch_startUpload):
        request = DummyRequest()
        request.matchdict['id'] = "26"
        request.POST['confirm'] = "true"
        request.user = createMockUser("username", "email", "password", 0)
        
        exp = createMockExperiment(26, 0, 0, 'preload')
        
        patch_getExperiment.return_value = exp
        
        result = upload_confirm_view(request)
        
        patch_getExperiment.assert_called_once_with(26, request.user, False)
        patch_startUpload.assert_called_once_with(exp)
        
        self.assertEqual(strings.experiment_upload_started_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_started_message % (request.application_url + "/account/experiments"), result['message'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual("true", result['confirm'])
        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_start_upload_view_should_get_confirmation(self, patch_getExperiment):
        request = DummyRequest()
        request.matchdict['id'] = "26"
        request.user = createMockUser("username", "email", "password", 0) 
        
        exp = createMockExperiment(26, 0, 0, 'preload')
        
        patch_getExperiment.return_value = exp
        
        result = upload_confirm_view(request)
        
        patch_getExperiment.assert_called_once_with(26, request.user, False)
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.experiment_upload_confirm_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_confirm_message, result['message'])
        self.assertEqual("false", result['confirm'])
        
        
        
class IntegrationTestUploadStatusView(IntegrationTestCase):
    def test_view_upload_already_started(self):
        from ptmscout.database import experiment
        self.bot.login()
        
        exp = experiment.getExperimentById(1, 0, False)
        exp.status = 'loading'
        exp.saveExperiment()
        
        result = self.ptmscoutapp.get("/upload/1", status=200)
        result.mustcontain(strings.experiment_upload_started_page_title)
    
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