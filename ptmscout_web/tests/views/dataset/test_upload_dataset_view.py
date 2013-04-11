from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.views.dataset import upload_dataset_view
from pyramid.testing import DummyRequest
from mock import patch, Mock
from ptmscout.config import strings, settings
from tests.views.mocking import createMockUser, createMockExperiment
from ptmscout.utils import forms
import os

class TestUploadDataset(UnitTestCase):
    
    @patch('ptmscout.views.dataset.upload_dataset_view.create_session')
    @patch('ptmscout.views.dataset.upload_dataset_view.create_experiment')
    @patch('ptmscout.utils.uploadutils.save_data_file')
    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.views.dataset.upload_dataset_view.create_schema')
    def test_view_should_handle_form_submission_redirect_on_success(self, patch_create_schema, patch_validator, patch_save, patch_create_exp, patch_create_session):
        request = DummyRequest()
        request.POST['datasetname'] = 'some name'
        request.POST['datasetfile'] = 'some file field'
        request.route_url = Mock()
        request.route_url.return_value = 'some_route' 
        
        request.user = createMockUser()
        exp = createMockExperiment(11, 0, 0)
        
        validator = patch_validator.return_value
        patch_create_schema.return_value = Mock(spec=forms.FormSchema())
        validator.validate.return_value = []

        patch_create_exp.return_value = 1001

        patch_save.return_value = "filename"
        patch_create_session.return_value = 111

        result = upload_dataset_view.upload_dataset_POST(request)
        self.assertEqual('some_route', result.location)

        request.route_url.assert_called_once_with('dataset_configure', id=111)

        patch_create_exp.assert_called_once_with('some name', "filename", request.user)
        patch_save.assert_called_once_with('some file field', settings.dataset_files_prefix)
        patch_create_session.assert_called_once_with(1001, "filename", request.user)
        
        validator.validate.assert_called_once_with()
        patch_create_schema.assert_called_once_with(request)
    
    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.views.dataset.upload_dataset_view.create_schema')
    def test_view_should_handle_form_submission_fail_on_check(self, patch_create_schema, patch_validator):
        request = DummyRequest()
        request.user = createMockUser()
        
        validator = patch_validator.return_value
        patch_create_schema.return_value = Mock(spec=forms.FormSchema())
        validator.validate.return_value = [strings.failure_reason_required_fields_cannot_be_empty]

        result = upload_dataset_view.upload_dataset_POST(request)

        validator.validate.assert_called_once_with()
        
        self.assertEqual(patch_create_schema.return_value, result['formrenderer'].schema)
        self.assertEqual([strings.failure_reason_required_fields_cannot_be_empty], result['errors'])
        self.assertEqual(strings.upload_annotations_page_title, result['pageTitle'])
        patch_create_schema.assert_called_once_with(request)
    
    @patch('ptmscout.views.dataset.upload_dataset_view.create_schema')
    def test_view_should_display(self, patch_create_schema):
        request = DummyRequest()
        request.user = createMockUser()
        
        patch_create_schema.return_value = Mock(spec=forms.FormSchema())
        
        result = upload_dataset_view.upload_dataset_GET(request)
        
        patch_create_schema.assert_called_once_with(request)
        
        self.assertEqual(patch_create_schema.return_value, result['formrenderer'].schema)
        self.assertEqual(strings.upload_annotations_page_title, result['pageTitle'])
        self.assertEqual([], result['errors'])
    

class IntegrationTestUploadDataset(IntegrationTestCase):

    def test_upload_file_should_start_session(self):
        result = self.ptmscoutapp.get("/dataset/upload", status=200)
        
        form = result.form
        filename = os.path.join(settings.ptmscout_path, "tests/behave/data/datasetLoad_correctDataset.txt")
        f = open(filename, 'rb')
        filecontents = f.read()
       
        form.set('datasetname', 'New dataset')
        form.set('datasetfile', (filename, filecontents))
        form.submit(status=302)
