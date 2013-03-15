from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.config import strings
from pyramid.testing import DummyRequest
from ptmscout.views.upload.upload_confirm import upload_confirm_view,\
    UploadAlreadyStarted, prepare_experiment
from tests.views.mocking import createMockExperiment, createMockUser,\
    createMockSession
from mock import patch
from ptmscout.database import upload
from ptmscout.utils import webutils

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

    def test_prepare_experiment_should_return_target_on_new(self):
        user = createMockUser()
        
        session = createMockSession(user)
        session.load_type = 'new'
        
        exp = createMockExperiment()
        exp.export = 0
        session.experiment_id = exp.id
        
        rval = prepare_experiment(session, exp, user)
        
        exp.saveExperiment.assert_called_once_with()
        
        self.assertEqual('in queue', exp.status)
        self.assertEqual(1, exp.export)
        self.assertEqual(exp, rval)
        
        self.assertEqual(exp.id, session.experiment_id)
        self.assertEqual('complete', session.stage)
        session.save.assert_called_once_with()
        
    def test_prepare_experiment_should_return_target_on_extend(self):
        user = createMockUser()
        
        session = createMockSession(user)
        session.load_type = 'extension'
        
        exp = createMockExperiment()
        exp.export = 0
        session.experiment_id = exp.id
        
        rval = prepare_experiment(session, exp, user)
        
        exp.saveExperiment.assert_called_once_with()
        
        self.assertEqual('in queue', exp.status)
        self.assertEqual(1, exp.export)
        self.assertEqual(exp, rval)
        
        self.assertEqual(exp.id, session.experiment_id)
        self.assertEqual('complete', session.stage)
        session.save.assert_called_once_with()

    @patch('ptmscout.database.experiment.getExperimentById')
    def test_prepare_experiment_should_copy_experiment_data_on_append(self, patch_getExperiment):
        user = createMockUser()
        
        session = createMockSession(user)
        session.load_type = 'append'
        
        exp = createMockExperiment()
        exp.export = 0
        session.experiment_id = exp.id
        
        target_exp = createMockExperiment()
        target_exp.export = 0
        session.parent_experiment = target_exp.id
        
        patch_getExperiment.return_value = target_exp
        
        rval = prepare_experiment(session, exp, user)
        
        target_exp.copyData.assert_called_once_with(exp)
        target_exp.saveExperiment.assert_called_once_with()
        exp.delete.assert_called_once_with()
        
        self.assertEqual('in queue', target_exp.status)
        self.assertEqual(0, exp.export)
        self.assertEqual(1, target_exp.export)
        self.assertEqual(target_exp, rval)
        
        self.assertEqual(target_exp.id, session.experiment_id)
        self.assertEqual('complete', session.stage)
        session.save.assert_called_once_with()
    
    @patch('ptmscout.database.modifications.deleteExperimentData')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_prepare_experiment_should_copy_experiment_data_and_delete_existing_on_reload(self, patch_getExperiment, patch_deleteData):
        user = createMockUser()
        
        session = createMockSession(user)
        session.load_type = 'reload'
        
        exp = createMockExperiment()
        exp.export = 0
        session.experiment_id = exp.id
        
        target_exp = createMockExperiment()
        target_exp.export = 0
        session.parent_experiment = target_exp.id
        
        patch_getExperiment.return_value = target_exp
        
        rval = prepare_experiment(session, exp, user)
        
        target_exp.copyData.assert_called_once_with(exp)
        patch_deleteData.assert_called_once_with(target_exp.id)
        
        target_exp.saveExperiment.assert_called_once_with()

        self.assertEqual('in queue', target_exp.status)
        self.assertEqual(0, exp.export)
        self.assertEqual(1, target_exp.export)
        self.assertEqual(target_exp, rval)     
        
        self.assertEqual(target_exp.id, session.experiment_id)
        self.assertEqual('complete', session.stage)
        session.save.assert_called_once_with()   

    @patch('ptmworker.data_import.start_import.apply_async')        
    @patch('ptmscout.views.upload.upload_confirm.prepare_experiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.upload.getSessionById')
    def test_start_upload_view_should_start_job_and_display_confirmation(self, patch_getSession, patch_getExperiment, patch_prepare, patch_startUpload):
        request = DummyRequest()
        request.POST['confirm'] = "true"
        request.POST['terms_of_use'] = "yes"
        request.matchdict['id'] = "102"
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=102, experiment_id=26, stage='confirm')
        exp = createMockExperiment(26, 0, None, 'configuration')
        target_exp = createMockExperiment(28, 0, None, 'preload')
        
        patch_getSession.return_value = session
        patch_getExperiment.return_value = exp
        exp.getUrl.return_value = "url"
        exp.getLongCitationString.return_value = "citation"
        
        patch_prepare.return_value = target_exp
        
        result = upload_confirm_view(request)
        
        patch_getSession.assert_called_once_with(102, request.user)
        patch_getExperiment.assert_called_once_with(26, request.user, False)

        patch_startUpload.assert_called_once_with((target_exp.id, session.id, request.user.email, request.application_url))
        
        exp_dict = webutils.object_to_dict(exp)
        exp_dict['url'] = "url"
        exp_dict['citation'] = "citation"
        self.assertEqual(exp_dict, result['experiment'])
        
        self.assertEqual(None, result['reason'])
        self.assertEqual(strings.experiment_upload_started_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_started_message % (request.application_url + "/account/experiments"), result['message'])
        self.assertEqual(102, result['session_id'])
        self.assertEqual(True, result['confirm'])
        
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.upload.getSessionById')
    def test_start_upload_view_should_fail_if_terms_not_accepted(self, patch_getSession, patch_getExperiment):
        request = DummyRequest()
        request.matchdict['id'] = "102"
        request.user = createMockUser() 
        request.POST['confirm'] = "true"
        
        session = createMockSession(request.user, sid=102, experiment_id=26, stage='confirm')
        exp = createMockExperiment(26, 0, None, 'in queue')
        
        patch_getSession.return_value=session
        patch_getExperiment.return_value = exp
        exp.getUrl.return_value = "url"
        exp.getLongCitationString.return_value = "citation"
        
        result = upload_confirm_view(request)
        
        patch_getSession.assert_called_once_with(102, request.user)
        patch_getExperiment.assert_called_once_with(26, request.user, False)
        
        exp_dict = webutils.object_to_dict(exp)
        exp_dict['url'] = "url"
        exp_dict['citation'] = "citation"
        self.assertEqual(exp_dict, result['experiment'])
        
        self.assertEqual(strings.failure_reason_terms_of_use_not_accepted, result['reason'])
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
        exp = createMockExperiment(26, 0, None, 'in queue')
        
        patch_getSession.return_value=session
        patch_getExperiment.return_value = exp
        exp.getUrl.return_value = "url"
        exp.getLongCitationString.return_value = "citation"
        
        result = upload_confirm_view(request)
        
        patch_getSession.assert_called_once_with(102, request.user)
        patch_getExperiment.assert_called_once_with(26, request.user, False)
        
        exp_dict = webutils.object_to_dict(exp)
        exp_dict['url'] = "url"
        exp_dict['citation'] = "citation"
        self.assertEqual(exp_dict, result['experiment'])
        
        self.assertEqual(None, result['reason'])
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
        session.change_name=''
        session.change_description=''
        session.stage='complete'
        session.save()
        
        result = self.ptmscoutapp.get("/upload/%d/confirm" % (session.id), status=200)
        result.mustcontain(strings.experiment_upload_started_page_title)
    
    @patch('ptmworker.data_import.start_import.apply_async')
    def test_view_integration(self, patch_startImport):
        from ptmscout.database import experiment
        self.bot.login()
        
        exp = experiment.getExperimentById(1, 0, False)
        exp.status = 'configuration'
        exp.saveExperiment()
        
        session = upload.Session()
        session.experiment_id = 1
        session.user_id = self.bot.user.id
        session.load_type='new'
        session.data_file='exp_file'
        session.parent_experiment=None
        session.change_name=''
        session.change_description=''
        session.stage='confirm'
        session.save()
        
        result = self.ptmscoutapp.get("/upload/%d/confirm" % session.id, status=200)
                
        result.mustcontain(strings.experiment_upload_confirm_message)
        result.form.set('terms_of_use', True)
        result = result.form.submit()
        
        result.mustcontain(strings.experiment_upload_started_message % ("http://localhost/account/experiments"))
        
        assert patch_startImport.called
