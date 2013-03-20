from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
from ptmscout.config import settings, strings
import os
from ptmscout.views.upload.upload_datafile import upload_data_file, create_session
from mock import patch, Mock
from tests.views.mocking import createMockUser, createMockExperiment,\
    createMockPermission
from pyramid.httpexceptions import HTTPFound
from ptmscout.utils import forms

class TestUploaddata_fileView(UnitTestCase):
    
    @patch('ptmscout.database.upload.Session')
    def test_create_session_should_create_new_session_and_return_session_id_when_extending(self, patch_session):
        request = DummyRequest()
        
        request.user = createMockUser()
        request.POST['load_type'] = 'extension'
        request.POST['parent_experiment'] = '24'
        request.POST['change_name'] = "Something Changed"
        request.POST['change_description'] = "some stuff"
        
        gen_session = patch_session.return_value
        gen_session.id = 10
        
        rval = create_session(request, 'exp_filename')
        
        gen_session.save.assert_called_once_with()
        
        self.assertEqual(gen_session.id, rval)
        
        self.assertEqual(request.user.id, gen_session.user_id)
        
        self.assertEqual('exp_filename', gen_session.data_file)
        self.assertEqual('extension', gen_session.load_type)
        self.assertEqual(24, gen_session.parent_experiment)
        self.assertEqual('Something Changed', gen_session.change_name)
        self.assertEqual('some stuff', gen_session.change_description)
        self.assertEqual('config', gen_session.stage)

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
        self.assertEqual('', gen_session.change_name)
        self.assertEqual('', gen_session.change_description)
        self.assertEqual('config', gen_session.stage)
    
    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.views.upload.upload_datafile.create_session')
    @patch('ptmscout.utils.uploadutils.save_data_file')
    @patch('ptmscout.views.upload.upload_datafile.create_schema')
    def test_view_should_save_data_file_create_session_and_forward_to_configuration(self, patch_create_schema, patch_save, patch_create_session, patch_validator):
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.POST['data_file'] = 'some file field'
        request.user = createMockUser()
        
        patch_create_schema.return_value = Mock(spec=forms.FormSchema())
        validator = patch_validator.return_value
        validator.validate.return_value = []
        
        patch_save.return_value = 'filename'
        patch_create_session.return_value = 1
        
        f = upload_data_file(request)
        assert isinstance(f, HTTPFound)
        
        self.assertEqual(request.application_url + '/upload/1/config', f.location)
        
        validator.validate.assert_called_once_with()
        patch_create_schema.assert_called_once_with(request, [])
        patch_save.assert_called_once_with('some file field', settings.experiment_files_prefix)
        patch_create_session.assert_called_once_with(request, 'filename')
    
    @patch('ptmscout.utils.forms.FormValidator')        
    @patch('ptmscout.views.upload.upload_datafile.create_schema')
    def test_view_should_handle_form_submission_fail_on_check(self, patch_create_schema, patch_validator):
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        validator = patch_validator.return_value
        patch_create_schema.return_value = Mock(spec=forms.FormSchema())
        validator.validate.return_value = [strings.failure_reason_required_fields_cannot_be_empty]

        result = upload_data_file(request)

        validator.validate.assert_called_once_with()
        
        self.assertEqual(patch_create_schema.return_value, result['formrenderer'].schema)
        self.assertEqual([strings.failure_reason_required_fields_cannot_be_empty], result['errors'])
        patch_create_schema.assert_called_once_with(request, [])
    
    @patch('ptmscout.views.upload.upload_datafile.create_schema')
    def test_view_should_get_users_experiments(self, patch_create_schema):
        request = DummyRequest()
        
        patch_create_schema.return_value = Mock(spec=forms.FormSchema())
        
        request.user = createMockUser()
        e1 = createMockExperiment(2, 0, 0)
        e1.ready.return_value = True
        e2 = createMockExperiment(3, 0, 0)
        e2.ready.return_value = True
        e3 = createMockExperiment(4, 0, 0)
        e3.ready.return_value = False
        p1 = createMockPermission(request.user, e1, 'owner')
        p2 = createMockPermission(request.user, e2, 'view')
        p3 = createMockPermission(request.user, e3, 'owner')
        request.user.permissions = [p1,p2,p3]
        
        result = upload_data_file(request)
        
        patch_create_schema.assert_called_once_with(request, [e1])
        
        self.assertEqual(patch_create_schema.return_value, result['formrenderer'].schema)
        self.assertEqual(strings.upload_page_title, result['pageTitle'])
        self.assertEqual([], result['errors'])
    

        
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
        form.set('change_name', "")
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
        form.set('change_name', "Changed")
        form.set('change_description', "This is a change")
        
        form.submit(status=302)
        
