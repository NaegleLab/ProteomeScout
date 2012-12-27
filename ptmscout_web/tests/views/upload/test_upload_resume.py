from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.database import upload
from tests.views.mocking import createMockSession, createMockUser
from pyramid.testing import DummyRequest
from pyramid.httpexceptions import HTTPFound
from ptmscout.views.upload.upload_resume import resume_upload_session
from mock import patch

class TestUploadResumeView(UnitTestCase):
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