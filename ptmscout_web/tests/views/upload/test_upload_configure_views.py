from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.database import upload
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockSession, createMockUser
from mock import patch
from ptmscout.config import strings, settings
import os
from pyramid.httpexceptions import HTTPFound
from ptmscout.utils.uploadutils import ColumnError, ErrorList
from ptmscout.views.upload.upload_configure import upload_config,\
    parse_user_input

class TestUploadConfigureView(UnitTestCase):

    def test_parse_user_input_should_accumulate_errors(self):
        request = DummyRequest()
        session = createMockSession(createMockUser())
        
        request.POST['submitted'] = "true"
        request.POST['override'] = "true"
        request.POST['units'] = ""
        
        request.POST['column_0_type'] = 'accession'
        
        request.POST['column_2_type'] = 'peptide'
        request.POST['column_3_type'] = 'modification'
        
        request.POST['column_4_type'] = 'data'
        request.POST['column_4_label'] = '0'
        
        request.POST['column_5_type'] = 'data'
        request.POST['column_5_label'] = ''
        
        request.POST['column_6_type'] = 'stddev'
        request.POST['column_6_label'] = '20'
        request.POST['column_7_type'] = 'stddev'
        request.POST['column_7_label'] = '20'
        
        request.POST['column_8_type'] = '   '

        
        column_defs, errors = parse_user_input(session, request)
        
        self.assertEqual([
                          strings.experiment_upload_error_column_type_not_defined % (2,),
                          strings.experiment_upload_error_data_column_empty_label % (6,),
                          strings.experiment_upload_error_data_column_label_duplicated % (8,),
                          strings.experiment_upload_error_column_type_not_defined % (9,),
                          strings.experiment_upload_error_standard_deviation_label_does_not_match_any_data_column % (7, "20"),
                          strings.experiment_upload_error_standard_deviation_label_does_not_match_any_data_column % (8, "20")
                          ], [ e.message for e in errors ])
        
        self.assertEqual({
                          'columns':{
                              0:{'type':'accession','label':''},
                              1:{'type':'', 'label':''},
                              2:{'type':'peptide','label':''},
                              3:{'type':'modification','label':''},
                              4:{'type':'data','label':'0'},
                              5:{'type':'data','label':''},
                              6:{'type':'stddev','label':'20'},
                              7:{'type':'stddev','label':'20'},
                              8:{'type':'', 'label':''}},
                          'units':''
                          }, column_defs)
        
    def test_parse_user_input_should_create_columns_for_all_form_elements(self):
        request = DummyRequest()
        session = createMockSession(createMockUser())
        
        request.POST['submitted'] = "true"
        request.POST['override'] = "true"
        request.POST['units'] = "time(min)"
        
        request.POST['column_0_type'] = 'accession'
        
        request.POST['column_1_type'] = 'hidden'
        
        request.POST['column_2_type'] = 'peptide'
        
        request.POST['column_3_type'] = 'modification'
        
        request.POST['column_4_type'] = 'run'
        
        request.POST['column_5_type'] = 'data'
        request.POST['column_5_label'] = '0'
        
        request.POST['column_6_type'] = 'data'
        request.POST['column_6_label'] = '20'
        
        request.POST['column_7_type'] = 'stddev'
        request.POST['column_7_label'] = '20'

        request.POST['column_8_type'] = 'none'
        request.POST['column_9_type'] = 'none'
        request.POST['column_9_label'] = 'some unused label'
        
        column_defs, errors = parse_user_input(session, request)
        
        self.maxDiff=None
        
        self.assertEqual([], errors)
        
        self.assertEqual({
                          'columns':{
                              0:{'type':'accession','label':''},
                              1:{'type':'hidden','label':''},
                              2:{'type':'peptide','label':''},
                              3:{'type':'modification','label':''},
                              4:{'type':'run','label':''},
                              5:{'type':'data','label':'0'},
                              6:{'type':'data','label':'20'},
                              7:{'type':'stddev','label':'20'},
                              8:{'type':'none','label':''},
                              9:{'type':'none','label':''}},
                          'units':'time(min)'
                          }, column_defs)
        
        self.assertEqual('accession', session.columns[0].type)
        self.assertEqual(0, session.columns[0].column_number)
        
        self.assertEqual('hidden', session.columns[1].type)
        self.assertEqual(1, session.columns[1].column_number)
        
        self.assertEqual('peptide', session.columns[2].type)
        self.assertEqual(2, session.columns[2].column_number)
        
        self.assertEqual('modification', session.columns[3].type)
        self.assertEqual(3, session.columns[3].column_number)
        
        self.assertEqual('run', session.columns[4].type)
        self.assertEqual(4, session.columns[4].column_number)
        
        self.assertEqual('data', session.columns[5].type)
        self.assertEqual('0', session.columns[5].label)
        self.assertEqual(5, session.columns[5].column_number)
        
        self.assertEqual('data', session.columns[6].type)
        self.assertEqual('20', session.columns[6].label)
        self.assertEqual(6, session.columns[6].column_number)
        
        self.assertEqual('stddev', session.columns[7].type)
        self.assertEqual('20', session.columns[7].label)
        self.assertEqual(7, session.columns[7].column_number)
        
        self.assertEqual('none', session.columns[8].type)
        self.assertEqual(8, session.columns[8].column_number)
        
        self.assertEqual('none', session.columns[9].type)
        self.assertEqual('', session.columns[9].label)
        self.assertEqual(9, session.columns[9].column_number)
        
        self.assertEqual('time(min)', session.units)
        
        
    @patch('ptmscout.views.upload.upload_configure.parse_user_input')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_stop_and_show_errors_disallow_force_if_error_in_form_input(self, patch_getSession, patch_parse):
        request = DummyRequest()
        request.matchdict['id'] = '234'
        request.POST['submitted'] = "true"
        request.POST['override'] = "true"
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=234)
        session.data_file = 'test/test_dataset_formatted.txt'
        session.load_type = 'new'
        session.stage = 'config'
        session.user_id = request.user.id
        patch_getSession.return_value = session
        
        column_vals = {"some":"defaults"}
        patch_parse.return_value = column_vals, [ColumnError("This is the error")]
        
        result = upload_config(request)

        patch_parse.assert_called_once_with(session, request)        
        patch_getSession.assert_called_once_with(234, request.user)
        
        expected_headers = open(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'test_dataset_formatted.txt'), 'r').readline().split("\t")
        expected_headers = [item for item in expected_headers if item.strip() != ""]
        
        self.assertEqual(False, result['allowoverride'])
        self.assertEqual(["This is the error"], result['error'])
        
        self.assertEqual(column_vals, result['data_definitions'])
        
        self.assertEqual(expected_headers, result['headers'])
        self.assertEqual(17, len(result['data_rows']))
        
        self.assertEqual(strings.experiment_upload_configure_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_configure_message, result['instruction'])


    @patch('ptmscout.views.upload.upload_configure.parse_user_input')
    @patch('ptmscout.utils.uploadutils.check_data_column_assignments')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_stop_and_show_errors_disallow_override_if_critical_error(self, patch_getSession, patch_check, patch_parse):
        request = DummyRequest()
        request.matchdict['id'] = '234'
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=234)
        session.data_file = 'test/test_dataset_formatted.txt'
        session.load_type = 'new'
        session.stage = 'config'
        session.user_id = request.user.id
        patch_getSession.return_value = session
        
        patch_check.side_effect = ErrorList([ColumnError("This is the error")], True)
        
        column_vals = {"some":"defaults"}
        patch_parse.return_value = column_vals, []
        
        result = upload_config(request)

        patch_parse.assert_called_once_with(session, request)        
        patch_getSession.assert_called_once_with(234, request.user)
        
        expected_headers = open(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'test_dataset_formatted.txt'), 'r').readline().split("\t")
        expected_headers = [item for item in expected_headers if item.strip() != ""]
        
        self.assertEqual(False, result['allowoverride'])
        self.assertEqual(["This is the error"], result['error'])
        
        self.assertEqual(column_vals, result['data_definitions'])
        
        self.assertEqual(expected_headers, result['headers'])
        self.assertEqual(17, len(result['data_rows']))
        
        self.assertEqual(strings.experiment_upload_configure_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_configure_message, result['instruction'])

    
    @patch('ptmscout.views.upload.upload_configure.parse_user_input')
    @patch('ptmscout.utils.uploadutils.check_data_column_assignments')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_stop_and_show_errors(self, patch_getSession, patch_check, patch_parse):
        request = DummyRequest()
        request.matchdict['id'] = '234'
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=234)
        session.data_file = 'test/test_dataset_formatted.txt'
        session.load_type = 'new'
        session.stage = 'config'
        session.user_id = request.user.id
        patch_getSession.return_value = session
        
        patch_check.side_effect = ErrorList([ColumnError("This is the error")], False)
        
        column_vals = {"some":"defaults"}
        patch_parse.return_value = column_vals, []
        
        result = upload_config(request)

        patch_parse.assert_called_once_with(session, request)        
        patch_getSession.assert_called_once_with(234, request.user)
        
        expected_headers = open(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'test_dataset_formatted.txt'), 'r').readline().split("\t")
        expected_headers = [item for item in expected_headers if item.strip() != ""]
        
        self.assertEqual(True, result['allowoverride'])
        self.assertEqual(["This is the error"], result['error'])
        
        self.assertEqual(column_vals, result['data_definitions'])
        
        self.assertEqual(expected_headers, result['headers'])
        self.assertEqual(17, len(result['data_rows']))
        
        self.assertEqual(strings.experiment_upload_configure_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_configure_message, result['instruction'])

    @patch('ptmscout.views.upload.upload_configure.parse_user_input')
    @patch('ptmscout.utils.uploadutils.check_data_column_assignments')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_configure_session_and_forward_to_metadata_even_with_errors_if_forced(self, patch_getSession, patch_check, patch_parse):
        request = DummyRequest()
        request.matchdict['id'] = '234'
        request.POST['submitted'] = "true"
        request.POST['override'] = "true"
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=234)
        session.data_file = 'test/test_dataset_formatted.txt'
        session.load_type = 'new'
        session.stage = 'config'
        session.user_id = request.user.id
        patch_getSession.return_value = session
        
        patch_parse.return_value = 'some_defs', []
        
        patch_check.side_effect = ErrorList([ColumnError("Ignored error")], False)
        
        f = upload_config(request)
        
        assert isinstance(f, HTTPFound)
        self.assertEqual(request.application_url + "/upload/%d/metadata" % (session.id), f.location)

        patch_parse.assert_called_once_with(session, request)
        patch_check.assert_called_once_with(session)
    
    @patch('ptmscout.views.upload.upload_configure.parse_user_input')
    @patch('ptmscout.utils.uploadutils.check_data_column_assignments')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_configure_session_and_forward_to_metadata(self, patch_getSession, patch_check, patch_parse):
        request = DummyRequest()
        request.matchdict['id'] = '234'
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=234)
        session.data_file = 'test/test_dataset_formatted.txt'
        session.load_type = 'new'
        session.stage = 'config'
        session.user_id = request.user.id
        patch_getSession.return_value = session
        patch_parse.return_value = 'some_defs', []
        
        f = upload_config(request)
        
        assert isinstance(f, HTTPFound)
        self.assertEqual(request.application_url + "/upload/%d/metadata" % (session.id), f.location)

        patch_parse.assert_called_once_with(session, request)
        patch_check.assert_called_once_with(session)
        session.save.assert_called_once_with()
    
    @patch('ptmscout.utils.uploadutils.assign_column_defaults')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_get_session_and_initialize_data(self, patch_getSession, patch_column_defaults):
        request = DummyRequest()
        request.matchdict['id'] = '234'
        request.user = createMockUser()
        
        session = createMockSession(request.user, sid=234)
        session.data_file = 'test/test_dataset_formatted.txt'
        session.load_type = 'new'
        session.stage = 'config'
        session.user_id = request.user.id
        patch_getSession.return_value = session
        
        def_column_vals = {"some":"defaults"}
        patch_column_defaults.return_value = def_column_vals
        
        result = upload_config(request)
        
        patch_getSession.assert_called_once_with(234, request.user)
        
        expected_headers = open(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'test_dataset_formatted.txt'), 'r').readline().split("\t")
        expected_headers = [item for item in expected_headers if item.strip() != ""]
        
        self.assertEqual(False, result['allowoverride'])
        self.assertEqual([], result['error'])
        
        self.assertEqual(def_column_vals, result['data_definitions'])
        self.assertEqual(expected_headers, result['headers'])
        self.assertEqual(17, len(result['data_rows']))
        
        self.assertEqual(strings.experiment_upload_configure_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_configure_message, result['instruction'])
        
        self.assertEqual(upload.SessionColumn.column_values, result['column_values'])
        
        
class IntegrationTestUploadConfigureView(IntegrationTestCase):
    def test_view_integration(self):
        self.bot.login()
        
        session = upload.Session()
        session.load_type='new'
        session.data_file='test/test_dataset_formatted.txt'
        session.stage='config'
        session.user_id=self.bot.user.id
        session.change_description=''
        session.save()
        
        self.ptmscoutapp.get("/upload/%d/config" % session.id, status=200)