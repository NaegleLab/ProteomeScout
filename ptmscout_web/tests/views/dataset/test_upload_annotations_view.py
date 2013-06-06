from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.views.dataset import upload_annotations_view
from tests.views.mocking import createMockExperiment, createMockUser
from pyramid.testing import DummyRequest
from ptmscout.utils import forms
from mock import patch, Mock
from ptmscout.config import strings, settings
import os

class TestAnnotateView(UnitTestCase):
    
    @patch('ptmscout.database.upload.Session')
    def test_create_session_should_save_annotation_resource_session(self, patch_session):
        session_obj = patch_session.return_value
        user = createMockUser()
        
        upload_annotations_view.create_session(12, 'name', 'filename', user)
        
        self.assertEqual('filename', session_obj.data_file)
        self.assertEqual(12, session_obj.experiment_id)
        self.assertEqual(user.id, session_obj.user_id)
        self.assertEqual('new', session_obj.load_type)
        self.assertEqual('annotations', session_obj.resource_type)
        self.assertEqual('config', session_obj.stage)
        self.assertEqual('name', session_obj.change_name)
        self.assertEqual('', session_obj.change_description)
        
        session_obj.save.assert_called_once_with()
    
    @patch('ptmscout.views.dataset.upload_annotations_view.create_session')
    @patch('ptmscout.utils.uploadutils.save_data_file')
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.utils.forms.FormValidator')        
    @patch('ptmscout.views.dataset.upload_annotations_view.create_schema')
    def test_view_should_handle_form_submission_redirect_on_success(self, patch_create_schema, patch_validator, patch_getExp, patch_save, patch_create_session):
        request = DummyRequest()
        request.POST['annotationname'] = 'some name'
        request.POST['annotationfile'] = 'some annotation file field'
        request.matchdict['id'] = '11'
        request.route_url = Mock()
        request.route_url.return_value = 'some_route' 
        
        request.user = createMockUser()
        exp = createMockExperiment(11, 0, 0)
        
        patch_getExp.return_value = exp
        validator = patch_validator.return_value
        mock_schema = Mock(spec=forms.FormSchema())
        patch_create_schema.return_value = mock_schema
        mock_schema.get_form_value.return_value = 'some name'

        validator.validate.return_value = []

        patch_save.return_value = "filename"
        patch_create_session.return_value = 111

        result = upload_annotations_view.upload_annotation_file_POST(request)
        self.assertEqual('some_route', result.location)

        request.route_url.assert_called_once_with('configure_annotations', id=11, sid=111)

        patch_save.assert_called_once_with('some annotation file field', settings.annotation_files_prefix)
        patch_create_session.assert_called_once_with(11, "some name", "filename", request.user)
        
        validator.validate.assert_called_once_with()
        patch_getExp.assert_called_once_with(exp.id, request.user)
        patch_create_schema.assert_called_once_with(request)
    
    
    
    
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.utils.forms.FormValidator')        
    @patch('ptmscout.views.dataset.upload_annotations_view.create_schema')
    def test_view_should_handle_form_submission_fail_on_check(self, patch_create_schema, patch_validator, patch_getExp):
        request = DummyRequest()
        request.matchdict['id'] = '11'
        request.user = createMockUser()
        exp = createMockExperiment(11, 0, 0)
        
        patch_getExp.return_value = exp
        validator = patch_validator.return_value
        patch_create_schema.return_value = Mock(spec=forms.FormSchema())
        validator.validate.return_value = [strings.failure_reason_required_fields_cannot_be_empty]

        result = upload_annotations_view.upload_annotation_file_POST(request)

        validator.validate.assert_called_once_with()
        
        patch_getExp.assert_called_once_with(exp.id, request.user)
        self.assertEqual(patch_create_schema.return_value, result['formrenderer'].schema)
        self.assertEqual([strings.failure_reason_required_fields_cannot_be_empty], result['errors'])
        self.assertEqual(strings.upload_annotations_page_title, result['pageTitle'])
        patch_create_schema.assert_called_once_with(request)
    
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.views.dataset.upload_annotations_view.create_schema')
    def test_view_should_get_users_experiments(self, patch_create_schema, patch_getExp):
        request = DummyRequest()
        request.matchdict['id'] = '11'
        request.user = createMockUser()
        exp = createMockExperiment(11, 0, 0)
        
        patch_create_schema.return_value = Mock(spec=forms.FormSchema())
        patch_getExp.return_value = exp
        
        result = upload_annotations_view.upload_annotation_file_GET(request)
        
        patch_create_schema.assert_called_once_with(request)
        
        patch_getExp.assert_called_once_with(exp.id, request.user)
        self.assertEqual(patch_create_schema.return_value, result['formrenderer'].schema)
        self.assertEqual(strings.upload_annotations_page_title, result['pageTitle'])
        self.assertEqual([], result['errors'])
        self.assertEqual(exp, result['experiment'])
    
        
class IntegrationTestAnnotateView(IntegrationTestCase):
    def test_upload_file_should_start_session(self):
        result = self.ptmscoutapp.get("/experiments/28/annotate", status=200)
        
        form = result.form
        filename = os.path.join(settings.ptmscout_path, "tests/behave/data/datasetLoad_correctDataset.txt")
        f = open(filename, 'rb')
        filecontents = f.read()
        
        form.set('annotationname', "some name")
        form.set('annotationfile', (filename, filecontents))
        form.submit(status=302)
