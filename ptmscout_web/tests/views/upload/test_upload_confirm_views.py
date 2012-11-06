from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.config import strings
from pyramid.testing import DummyRequest
from ptmscout.views.upload.upload_confirm import upload_confirm_view,\
    UploadAlreadyStarted
from tests.views.mocking import createMockExperiment, createMockUser,\
    createMockSession
from mock import patch
from ptmscout.database import upload

class TestUploadStatusView(UnitTestCase):
    
    @patch('ptmscout.database.upload.getSessionById')
    def test_start_upload_view_should_redirect_if_already_started(self, patch_getSession):
        request = DummyRequest()
        request.matchdict['id'] = "102"
        request.user = createMockUser() 
        
        session = createMockSession(request.user, sid=102, experiment_id=26, stage='complete')
        patch_getSession.return_value = session

        try:
            upload_confirm_view(request)
        except UploadAlreadyStarted:
            pass
        else:
            self.fail("Expected exception UploadAlreadyStarted")
    
    @patch('ptmworker.tasks.start_import.apply_async')        
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.upload.getSessionById')
    def test_start_upload_view_should_start_job_and_display_confirmation(self, patch_getSession, patch_getExperiment, patch_startUpload):
        request = DummyRequest()
        request.POST['confirm'] = "true"
        request.POST['terms_of_use'] = "yes"
        request.matchdict['id'] = "102"
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=102, experiment_id=26, stage='confirm')
        exp = createMockExperiment(26, 0, None, 'preload')
        
        patch_getSession.return_value = session
        patch_getExperiment.return_value = exp
        
        result = upload_confirm_view(request)
        
        patch_getSession.assert_called_once_with(102, request.user)
        patch_getExperiment.assert_called_once_with(26, request.user, False)
        session.save.assert_called_once_with()
        patch_startUpload.assert_called_once_with((request, exp, session, request.user))
        
        self.assertEqual('complete', session.stage)
        self.assertEqual(None, result['reason'])
        self.assertEqual(strings.experiment_upload_started_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_started_message % (request.application_url + "/account/experiments"), result['message'])
        self.assertEqual(102, result['session_id'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(True, result['confirm'])
        
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.upload.getSessionById')
    def test_start_upload_view_should_fail_if_terms_not_accepted(self, patch_getSession, patch_getExperiment):
        request = DummyRequest()
        request.matchdict['id'] = "102"
        request.user = createMockUser() 
        request.POST['confirm'] = "true"
        
        session = createMockSession(request.user, sid=102, experiment_id=26, stage='confirm')
        exp = createMockExperiment(26, 0, None, 'preload')
        
        patch_getSession.return_value=session
        patch_getExperiment.return_value = exp
        
        result = upload_confirm_view(request)
        
        patch_getSession.assert_called_once_with(102, request.user)
        patch_getExperiment.assert_called_once_with(26, request.user, False)
        
        self.assertEqual(strings.failure_reason_terms_of_use_not_accepted, result['reason'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.experiment_upload_confirm_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_confirm_message, result['message'])
        self.assertEqual(102, result['session_id'])
        self.assertEqual(False, result['confirm'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.upload.getSessionById')
    def test_start_upload_view_should_get_confirmation(self, patch_getSession, patch_getExperiment):
        request = DummyRequest()
        request.matchdict['id'] = "102"
        request.user = createMockUser() 
        
        session = createMockSession(request.user, sid=102, experiment_id=26, stage='confirm')
        exp = createMockExperiment(26, 0, None, 'preload')
        
        patch_getSession.return_value=session
        patch_getExperiment.return_value = exp
        
        result = upload_confirm_view(request)
        
        patch_getSession.assert_called_once_with(102, request.user)
        patch_getExperiment.assert_called_once_with(26, request.user, False)
        self.assertEqual(None, result['reason'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.experiment_upload_confirm_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_confirm_message, result['message'])
        self.assertEqual(102, result['session_id'])
        self.assertEqual(False, result['confirm'])
        
        
        
class IntegrationTestUploadStatusView(IntegrationTestCase):
    def test_view_upload_already_started(self):
        from ptmscout.database import experiment
        self.bot.login()
        
        exp = experiment.getExperimentById(1, 0, False)
        exp.status = 'loading'
        exp.saveExperiment()
        
        session = upload.Session()
        session.experiment_id = 1
        session.user_id = self.bot.user.id
        session.load_type='new'
        session.data_file='exp_file'
        session.parent_experiment=None
        session.change_description=''
        session.stage='complete'
        session.save()
        
        result = self.ptmscoutapp.get("/upload/%d/confirm" % (session.id), status=200)
        result.mustcontain(strings.experiment_upload_started_page_title)
    
    @patch('ptmworker.tasks.start_import.apply_async')
    def test_view_integration(self, patch_startImport):
        from ptmscout.database import experiment
        self.bot.login()
        
        exp = experiment.getExperimentById(1, 0, False)
        exp.status = 'preload'
        exp.saveExperiment()
        
        session = upload.Session()
        session.experiment_id = 1
        session.user_id = self.bot.user.id
        session.load_type='new'
        session.data_file='exp_file'
        session.parent_experiment=None
        session.change_description=''
        session.stage='confirm'
        session.save()
        
        result = self.ptmscoutapp.get("/upload/%d/confirm" % session.id, status=200)
                
        result.mustcontain(strings.experiment_upload_confirm_message)
        result.form.set('terms_of_use', True)
        result = result.form.submit()
        
        result.mustcontain(strings.experiment_upload_started_message % ("http://localhost/account/experiments"))
        
        assert patch_startImport.called