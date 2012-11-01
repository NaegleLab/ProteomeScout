from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
from ptmscout.config import settings, strings
import os
from ptmscout.views.upload.upload_datafile import save_data_file,upload_data_file,\
    check_required_fields, create_session
from mock import patch
import cgi
from tests.views.mocking import createMockUser, createMockExperiment,\
    createMockPermission
from pyramid.httpexceptions import HTTPFound

class TestUploaddata_fileView(UnitTestCase):
    
    @patch('ptmscout.database.upload.Session')
    def test_create_session_should_create_new_session_and_return_session_id(self, patch_session):
        request = DummyRequest()
        
        request.user = createMockUser()
        request.POST['load_type'] = 'new'
        request.POST['parent_experiment'] = '24'
        request.POST['change_description'] = "some stuff"
        
        gen_session = patch_session.return_value
        gen_session.id = 10
        
        rval = create_session(request, 'exp_filename')
        
        gen_session.save.assert_called_once_with()
        
        self.assertEqual(gen_session.id, rval)
        
        self.assertEqual(request.user.id, gen_session.user_id)
        
        self.assertEqual('exp_filename', gen_session.data_file)
        self.assertEqual('new', gen_session.load_type)
        self.assertEqual(None, gen_session.parent_experiment)
        self.assertEqual('', gen_session.change_description)
        self.assertEqual('config', gen_session.stage)
        
        
    @patch('ptmscout.views.upload.upload_datafile.verify_data_file')
    def test_save_data_file_should_save_experiment_data_file_and_pass(self, patch_verify):
        request = DummyRequest()
        request.POST['data_file'] = cgi.FieldStorage()
        
        cwd = settings.ptmscout_path
        os.chdir(cwd)
        
        filename = "datasetLoad_correctDataset.txt"
        request.POST['data_file'].filename = filename
        request.POST['data_file'].file = open(os.sep.join(["tests","behave","data",filename]), 'rb')
        
        patch_verify.return_value = None
        
        status, exp_file = save_data_file(request)
        
        os.chdir(os.sep.join(["tests", "behave", "data"]))
        os.system("mac2unix -q -n %s tmp_test.conv" % (filename))
        os.system("dos2unix -q tmp_test.conv")

        tf = open("tmp_test.conv", 'rb')
        test_file = tf.read()
        tf.close()
        
        os.remove("tmp_test.conv")
        os.chdir(cwd)
        
        os.chdir(settings.experiment_data_file_path)
        tf = open(exp_file, 'rb')
        result_file = tf.read()
        tf.close()
        
        os.remove(exp_file)
        
        patch_verify.assert_called_once_with(exp_file)
        
        self.assertEqual(True, status)
        self.assertEqual(test_file, result_file)
    
    
    def test_check_required_fields_should_return_false_if_extension_experiment_and_no_change(self):
        request = DummyRequest()
        
        data_file_field = cgi.FieldStorage()
        request.POST['load_type'] = "extension"
        request.POST['data_file'] = data_file_field 
        request.POST['parent_experiment'] = "24"
        request.POST['change_description'] = "   "
        
        status, error, field_dict = check_required_fields(request, [createMockExperiment(eid=24)])

        self.assertEqual(False, status)
        self.assertEqual(strings.failure_reason_required_fields_cannot_be_empty % "Change Description", error)
        self.assertEqual({'load_type':"extension", 'data_file': data_file_field, 'parent_experiment':"24", 'change_description':""}, field_dict)

    def test_check_required_fields_should_return_false_if_derivative_experiment_and_bad_parent_selected(self):
        request = DummyRequest()
        
        data_file_field = cgi.FieldStorage()
        request.POST['load_type'] = "append"
        request.POST['data_file'] = data_file_field 
        request.POST['parent_experiment'] = "20"
        
        status, error, field_dict = check_required_fields(request, [createMockExperiment(eid=24)])

        self.assertEqual(False, status)
        self.assertEqual(strings.failure_reason_field_value_not_valid % "Parent Experiment", error)
        self.assertEqual({'load_type':"append",'data_file': data_file_field, 'parent_experiment':"20"}, field_dict)

    
    def test_check_required_fields_should_return_false_if_derivative_experiment_and_no_parent_selected(self):
        request = DummyRequest()
        
        data_file_field = cgi.FieldStorage()
        request.POST['load_type'] = "reload"
        request.POST['data_file'] = data_file_field 
        request.POST['parent_experiment'] = ""
        request.POST['change_description'] = "Some changes"
        
        status, error, field_dict = check_required_fields(request, [])

        self.assertEqual(False, status)
        self.assertEqual(strings.failure_reason_required_fields_cannot_be_empty % "Parent Experiment", error)
        self.assertEqual({'load_type':"reload", 'data_file': data_file_field, 'parent_experiment':"", 'change_description':"Some changes"}, field_dict)

    
    @patch('ptmscout.views.upload.upload_datafile.save_data_file')
    @patch('ptmscout.views.upload.upload_datafile.check_required_fields')
    def test_view_should_handle_form_submission_fail_on_data_file_failure(self, patch_check, patch_verify):
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        field_dict = {'name':"something"}
        patch_check.return_value = (True, None, field_dict)
        patch_verify.return_value = (False, strings.failure_reason_experiment_file_not_enough_columns)
        
        result = upload_data_file(request)
        
        self.assertEqual(strings.failure_reason_experiment_file_not_enough_columns, result['reason'])
        self.assertEqual(field_dict, result['formfields'])
        patch_check.assert_called_once_with(request, [])
        patch_verify.assert_called_once_with(request)
    
    @patch('ptmscout.views.upload.upload_datafile.create_session')
    @patch('ptmscout.views.upload.upload_datafile.save_data_file')
    @patch('ptmscout.views.upload.upload_datafile.check_required_fields')
    def test_view_should_save_data_file_create_session_and_forward_to_configuration(self, patch_check, patch_save, patch_create_session):
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        field_dict = {'name':"something"}
        patch_check.return_value = (True, None, field_dict)
        patch_save.return_value = (True, 'filename')
        patch_create_session.return_value = 1
        
        f = upload_data_file(request)
        assert isinstance(f, HTTPFound)
        
        self.assertEqual(request.application_url + '/upload/1/config', f.location)
        patch_check.assert_called_once_with(request, [])
        patch_save.assert_called_once_with(request)
        patch_create_session.assert_called_once_with(request, 'filename')
        
    @patch('ptmscout.views.upload.upload_datafile.check_required_fields')
    def test_view_should_handle_form_submission_fail_on_check(self, patch_check):
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        field_dict = {'name':"something"}
        patch_check.return_value = (False, strings.failure_reason_required_fields_cannot_be_empty, field_dict)

        result = upload_data_file(request)

        self.assertEqual(strings.failure_reason_required_fields_cannot_be_empty, result['reason'])
        self.assertEqual(field_dict, result['formfields'])
        patch_check.assert_called_once_with(request, [])
    
    def test_view_should_get_users_experiments(self):
        request = DummyRequest()
        
        request.user = createMockUser()
        e1 = createMockExperiment(2, 0, 0)
        e2 = createMockExperiment(3, 0, 0)
        p1 = createMockPermission(request.user, e1, 'owner')
        p2 = createMockPermission(request.user, e2, 'view')
        request.user.permissions = [p1,p2]
        
        result = upload_data_file(request)
        
        expected_experiments = [{'method_calls':[], 'status':'loaded', 'URL':e1.URL, 'parent_id': 0, 'id': 2, 'name': 'Experiment Name2', 'public': 0}]
        
        del result['user_experiments'][0]['date']
        
        self.assertEqual(strings.upload_page_title, result['pageTitle'])
        
        self.assertEqual(expected_experiments, result['user_experiments'])
        self.assertEqual({'load_type':"", 'parent_experiment':"", 'change_description':""}, result['formfields'])
        self.assertEqual(None, result['reason'])
    

        
class IntegrationTestUploadView(IntegrationTestCase):
    
    def test_forbidden_should_invoke_on_unauthorized_access(self):
        self.bot.logout()
        response = self.ptmscoutapp.get('/upload')
        response.mustcontain("forbidden")
    
    def test_view_integration_when_missing_fields(self):
        result = self.ptmscoutapp.get("/upload", status=200)
        
        form = result.form
        
        form.set('load_type', "new")
        
        result = form.submit(status=200)
        result.mustcontain(strings.failure_reason_required_fields_cannot_be_empty % "Input Data File")
    
    def test_view_integration(self):
        self.bot.acquire_experiments([26])
        result = self.ptmscoutapp.get("/upload", status=200)
        result.mustcontain(self.bot.user.permissions[0].experiment.name)
        
        form = result.form
        filename = os.path.join(settings.ptmscout_path, "tests/behave/data/datasetLoad_correctDataset.txt")
        f = open(filename, 'rb')
        filecontents = f.read()
        
        form.set('data_file', (filename, filecontents))
        form.set('load_type', "append")
        form.set('parent_experiment', "26")
        form.set('change_description', "")
        
        form.submit(status=302)
        
    def test_view_integration_when_extending(self):
        self.bot.acquire_experiments([26])
        result = self.ptmscoutapp.get("/upload", status=200)
        result.mustcontain(self.bot.user.permissions[0].experiment.name)
        
        form = result.form
        filename = os.path.join(settings.ptmscout_path, "tests/behave/data/datasetLoad_correctDataset.txt")
        f = open(filename, 'rb')
        filecontents = f.read()
        
        form.set('data_file', (filename, filecontents))
        form.set('load_type', "extension")
        form.set('parent_experiment', "26")
        form.set('change_description', "This is a change")
        
        form.submit(status=302)
        