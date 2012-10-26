from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.database import upload
from pyramid.testing import DummyRequest
from ptmscout.views.upload.upload_configure import upload_config, ColumnError
from tests.views.mocking import createMockSession, createMockUser
from mock import patch
from ptmscout.config import strings, settings
import os
from pyramid.httpexceptions import HTTPFound

class TestUploadConfigureView(UnitTestCase):

    @patch('ptmscout.views.upload.upload_configure.assign_column_defaults')
    @patch('ptmscout.views.upload.upload_configure.check_data_column_assignments')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_stop_and_show_errors(self, patch_getSession, patch_check, patch_column_defaults):
        request = DummyRequest()
        request.matchdict['id'] = '234'
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=234)
        session.data_file = 'test/test_dataset.txt'
        session.load_type = 'new'
        session.stage = 'config'
        session.user_id = request.user.id
        patch_getSession.return_value = session
        
        patch_check.side_effect = ColumnError("This is the error")
        result = upload_config(request)
        
        def_column_vals = {"some":"defaults"}
        patch_column_defaults.return_value = def_column_vals
        
        patch_getSession.assert_called_once_with(234, request.user)
        
        expected_headers = open(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'test_dataset.txt'), 'r').readline().split("\t")
        
        self.assertEqual("This is the error", result['error'])
        
        self.assertEqual(def_column_vals, result['columns'])
        self.assertEqual(expected_headers, result['headers'])
        self.assertEqual(18, len(result['data_rows']))
        
        self.assertEqual(strings.experiment_upload_configure_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_configure_message, result['instruction'])

    @patch('ptmscout.views.upload.upload_configure.check_data_column_assignments')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_configure_session_and_forward_to_metadata_even_with_errors_if_forced(self, patch_getSession, patch_check):
        request = DummyRequest()
        request.matchdict['id'] = '234'
        request.POST['submitted'] = "true"
        request.POST['override'] = "true"
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=234)
        session.data_file = 'test/test_dataset.txt'
        session.load_type = 'new'
        session.stage = 'config'
        session.user_id = request.user.id
        patch_getSession.return_value = session
        
        patch_check.side_effect = ColumnError("")
        
        f = upload_config(request)
        
        assert isinstance(f, HTTPFound)
        self.assertEqual(request.application_url + "/upload/%d/metadata" % (), f.location)

        patch_check.assert_called_once_with(session)
    
    @patch('ptmscout.views.upload.upload_configure.check_data_column_assignments')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_configure_session_and_forward_to_metadata(self, patch_getSession, patch_check):
        request = DummyRequest()
        request.matchdict['id'] = '234'
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=234)
        session.data_file = 'test/test_dataset.txt'
        session.load_type = 'new'
        session.stage = 'config'
        session.user_id = request.user.id
        patch_getSession.return_value = session
        
        f = upload_config(request)
        
        assert isinstance(f, HTTPFound)
        self.assertEqual(request.application_url + "/upload/%d/metadata" % (), f.location)

        patch_check.assert_called_once_with(session)
    
    @patch('ptmscout.views.upload.upload_configure.assign_column_defaults')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_get_session_and_initialize_data(self, patch_getSession, patch_column_defaults):
        request = DummyRequest()
        request.matchdict['id'] = '234'
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=234)
        session.data_file = 'test/test_dataset.txt'
        session.load_type = 'new'
        session.stage = 'config'
        session.user_id = request.user.id
        patch_getSession.return_value = session
        
        def_column_vals = {"some":"defaults"}
        patch_column_defaults.return_value = def_column_vals
        
        result = upload_config(request)
        
        patch_getSession.assert_called_once_with(234, request.user)
        
        expected_headers = open(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'test_dataset.txt'), 'r').readline().split("\t")
        
        self.assertEqual(None, result['error'])
        
        self.assertEqual(def_column_vals, result['columns'])
        self.assertEqual(expected_headers, result['headers'])
        self.assertEqual(18, len(result['data_rows']))
        
        self.assertEqual(strings.experiment_upload_configure_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_configure_message, result['instruction'])
        
        
class IntegrationTestUploadConfigureView(IntegrationTestCase):
    def test_view_integration(self):
        self.bot.login()
        
        session = upload.Session()
        session.load_type='new'
        session.data_file='test/test_dataset.txt'
        session.stage='config'
        session.user_id=self.bot.user.id
        session.change_description=''
        session.save()
        
        self.ptmscoutapp.get("/upload/%d/config" % session.id, status=200)