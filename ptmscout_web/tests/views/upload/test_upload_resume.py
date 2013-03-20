from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.database import upload
from tests.views.mocking import createMockSession, createMockUser, createMockExperiment
from pyramid.testing import DummyRequest
from pyramid.httpexceptions import HTTPFound
from ptmscout.views.upload.upload_resume import resume_upload_session, retry_failed_upload
from ptmscout.views.upload import upload_confirm
from mock import patch
from ptmscout.utils import webutils
from ptmscout.config import strings

class TestUploadResumeView(UnitTestCase):

    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.upload.getSessionById')
    def test_retry_upload_view_should_redirect_if_already_started(self, patch_getSession, patch_getExperiment):
        request = DummyRequest()
        request.matchdict['id'] = "102"
        request.user = createMockUser() 
        
        session = createMockSession(request.user, sid=102, experiment_id=26, stage='complete')
        patch_getSession.return_value = session

        exp = createMockExperiment(26, 0, None, 'loading')
        patch_getExperiment.return_value = exp

        try:
            retry_failed_upload(request)
        except upload_confirm.UploadAlreadyStarted:
            pass
        else:
            self.fail("Expected exception UploadAlreadyStarted")


    @patch('ptmworker.data_import.start_import.apply_async')        
    @patch('ptmscout.views.upload.upload_resume.prepare_experiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.upload.getSessionById')
    def test_retry_upload_view_should_start_job_and_display_confirmation(self, patch_getSession, patch_getExperiment, patch_prepare, patch_startUpload):
        request = DummyRequest()
        request.POST['confirm'] = "true"
        request.POST['terms_of_use'] = "yes"
        request.matchdict['id'] = "102"
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=102, experiment_id=26, stage='confirm')
        exp = createMockExperiment(26, 0, None, 'error')
        
        patch_getSession.return_value = session
        patch_getExperiment.return_value = exp
        exp.getUrl.return_value = "url"
        exp.getLongCitationString.return_value = "citation"
        
        result = retry_failed_upload(request)
        
        patch_getSession.assert_called_once_with(102, request.user)
        patch_getExperiment.assert_called_once_with(26, request.user, False)

        patch_startUpload.assert_called_once_with((exp.id, session.id, exp.job.id))
        
        exp_dict = webutils.object_to_dict(exp)
        exp_dict['url'] = "url"
        exp_dict['citation'] = "citation"
        self.assertEqual(exp_dict, result['experiment'])
        
        self.assertEqual(None, result['reason'])
        self.assertEqual(strings.experiment_upload_started_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_started_message % (request.application_url + "/account/experiments"), result['message'])
        self.assertEqual(102, result['session_id'])
        self.assertEqual(True, result['confirm'])
        


    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_redirect_to_config(self, patch_getSession):
        user = createMockUser()
        session = createMockSession(user)
        session.stage = 'config'
        
        request = DummyRequest()
        request.matchdict['id'] = str(session.id)
        request.user = user
        
        patch_getSession.return_value = session
        
        f = resume_upload_session(request)
        
        self.assertTrue(isinstance(f, HTTPFound))
        self.assertEqual(request.application_url + "/upload/%d/config" % (session.id), f.location)
                
        patch_getSession.assert_called_once_with(session.id, user)
        
    
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_redirect_to_metadata(self, patch_getSession):
        user = createMockUser()
        session = createMockSession(user)
        session.stage = 'metadata'
        
        request = DummyRequest()
        request.matchdict['id'] = str(session.id)
        request.user = user
        
        patch_getSession.return_value = session
        
        f = resume_upload_session(request)
        
        self.assertTrue(isinstance(f, HTTPFound))
        self.assertEqual(request.application_url + "/upload/%d/metadata" % (session.id), f.location)
                
        patch_getSession.assert_called_once_with(session.id, user)
        
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_redirect_to_conditions(self, patch_getSession):
        user = createMockUser()
        session = createMockSession(user)
        session.stage = 'condition'
        
        request = DummyRequest()
        request.matchdict['id'] = str(session.id)
        request.user = user
        
        patch_getSession.return_value = session
        
        f = resume_upload_session(request)
        
        self.assertTrue(isinstance(f, HTTPFound))
        self.assertEqual(request.application_url + "/upload/%d/conditions" % (session.id), f.location)
                
        patch_getSession.assert_called_once_with(session.id, user)
        
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_redirect_to_confirm(self, patch_getSession):
        user = createMockUser()
        session = createMockSession(user)
        session.stage = 'confirm'
        
        request = DummyRequest()
        request.matchdict['id'] = str(session.id)
        request.user = user
        
        patch_getSession.return_value = session
        
        f = resume_upload_session(request)
        
        self.assertTrue(isinstance(f, HTTPFound))
        self.assertEqual(request.application_url + "/upload/%d/confirm" % (session.id), f.location)
                
        patch_getSession.assert_called_once_with(session.id, user)    
        
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_redirect_to_confirm_if_complete(self, patch_getSession):
        user = createMockUser()
        session = createMockSession(user)
        session.stage = 'complete'
        
        request = DummyRequest()
        request.matchdict['id'] = str(session.id)
        request.user = user
        
        patch_getSession.return_value = session
        
        f = resume_upload_session(request)
        
        self.assertTrue(isinstance(f, HTTPFound))
        self.assertEqual(request.application_url + "/upload/%d/confirm" % (session.id), f.location)
                
        patch_getSession.assert_called_once_with(session.id, user)
        
class IntegrationTestUploadResumeView(IntegrationTestCase):
    def test_view_should_forbid(self):
        session = upload.Session()
        session.experiment_id = 26
        session.data_file = ''
        session.change_name = ''
        session.change_description = ''
        session.load_type = 'new'
        session.stage = 'confirm'
        session.units = ''
        session.user_id = self.bot.user.id
        session.save()
        
        self.bot.logout()
        
        result = self.ptmscoutapp.get("/upload/%d" % (session.id), status=200)
        result.mustcontain("Forbidden")
    
    def test_view_integration(self):
        session = upload.Session()
        session.experiment_id = 26
        session.data_file = ''
        session.change_name = ''
        session.change_description = ''
        session.load_type = 'new'
        session.stage = 'confirm'
        session.units = ''
        session.user_id = self.bot.user.id
        session.save()
        
        self.bot.login()
        
        result = self.ptmscoutapp.get("/upload/%d" % (session.id), status=302)
        nresult = result.follow(status=200)
        
        nresult.mustcontain("confirm")
