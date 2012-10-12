from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
from ptmscout.views.experiment.upload import user_upload, check_required_fields
from ptmscout.config import strings
from tests.views.mocking import createMockUser, createMockPermission,\
    createMockExperiment
import base64
import json
from mock import patch

class TestUploadView(UnitTestCase):

    @patch('ptmscout.views.experiment.upload.get_enum_fields')
    @patch('ptmscout.views.experiment.upload.get_numeric_fields')
    @patch('ptmscout.views.experiment.upload.get_required_fields')
    def test_check_required_fields_should_fail_on_enum_conflict(self, patch_required, patch_numeric, patch_enum):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "yes"
        request.POST['parent_experiment'] = "100"
        request.POST['publication_year'] = "2008"
        request.POST['terms_of_use'] = "yes"
        
        patch_required.return_value = ['published']
        patch_numeric.return_value = ['publication_year', 'pmid']
        patch_enum.return_value = {'parent_experiment': set(["24","26"]), 'published':set(["yes","no"])}
        
        exp_ids = [24, 26]
        status, error, field_dict = check_required_fields(request, exp_ids)
        
        patch_required.assert_called_once_with(request)
        patch_numeric.assert_called_once_with(request)
        patch_enum.assert_called_once_with(request, exp_ids)
        
        self.assertEqual(False, status)
        self.assertEqual(strings.failure_reason_field_value_not_valid % "Parent Experiment", error)
        self.assertEqual({'pmid':"1234", 'URL':"http://somestuff.com", 'published':"yes", 'parent_experiment':"100", 'publication_year':"2008"}, field_dict)
        

    @patch('ptmscout.views.experiment.upload.get_numeric_fields')
    @patch('ptmscout.views.experiment.upload.get_required_fields')
    def test_check_required_fields_should_fail_on_numeric_field_format_incorrect(self, patch_required, patch_numeric):
        request = DummyRequest()
        request.POST['pmid'] = "na1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "yes"
        request.POST['parent_experiment'] = "24"
        request.POST['publication_year'] = "2008"
        request.POST['terms_of_use'] = "yes"
        
        patch_required.return_value = ['published']
        patch_numeric.return_value = ['publication_year', 'pmid']
        
        exp_ids = [24, 26]
        status, error, field_dict = check_required_fields(request, exp_ids)
        
        patch_required.assert_called_once_with(request)
        patch_numeric.assert_called_once_with(request)
        
        self.assertEqual(False, status)
        self.assertEqual(strings.failure_reason_field_must_be_numeric % "PubMed ID", error)
        self.assertEqual({'pmid':"na1234", 'URL':"http://somestuff.com", 'published':"yes", 'parent_experiment':"24", 'publication_year':"2008"}, field_dict)
        

    @patch('ptmscout.views.experiment.upload.get_required_fields')
    def test_check_required_fields_should_return_false_if_empty(self, patch_required):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "   "
        request.POST['terms_of_use'] = "yes"
        
        patch_required.return_value = ['published']
        
        status, error, field_dict = check_required_fields(request, [])
        
        patch_required.assert_called_once_with(request)
        
        self.assertEqual(False, status)
        self.assertEqual(strings.failure_reason_required_fields_cannot_be_empty, error)
        self.assertEqual({'pmid':"1234", 'URL':"http://somestuff.com", 'published':""}, field_dict)

    def test_check_required_fields_should_return_false_if_terms_not_accepted(self):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "   "
        
        status, error, field_dict = check_required_fields(request, [])
        
        self.assertEqual(False, status)
        self.assertEqual(strings.failure_reason_terms_of_use_not_accepted, error)
        self.assertEqual({'pmid':"1234", 'URL':"http://somestuff.com", 'published':""}, field_dict)


    @patch('ptmscout.views.experiment.upload.get_numeric_fields')
    @patch('ptmscout.views.experiment.upload.get_enum_fields')
    @patch('ptmscout.views.experiment.upload.get_required_fields')
    def test_check_required_fields_should_return_true_on_success(self, patch_required, patch_enums, patch_numeric):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "yes"
        request.POST['parent_experiment'] = "24"
        request.POST['publication_year'] = "2008"
        request.POST['terms_of_use'] = "yes"
        
        patch_required.return_value = ['published']
        patch_numeric.return_value = ['publication_year']
        patch_enums.return_value = {'parent_experiment':set(["24", "26"])}
        
        exp_ids = [24, 26]
        status, error, field_dict = check_required_fields(request, exp_ids)
        
        patch_required.assert_called_once_with(request)
        patch_numeric.assert_called_once_with(request)
        patch_enums.assert_called_once_with(request, exp_ids)
        
        self.assertEqual(True, status)
        self.assertEqual(None, error)
        self.assertEqual({'pmid':"1234", 'URL':"http://somestuff.com", 'published':"yes", 'parent_experiment':"24", 'publication_year':"2008"}, field_dict)
        

    @patch('ptmscout.views.experiment.upload.create_experiment_and_start_upload')
    @patch('ptmscout.views.experiment.upload.verify_datafile')
    @patch('ptmscout.views.experiment.upload.check_required_fields')
    def test_view_should_start_upload_on_successful_verification(self, patch_check, patch_verify, patch_start):
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.user = createMockUser("username", "email", "password", 1)
        
        field_dict = {'name':"something"}
        patch_check.return_value = (True, None, field_dict)
        tmp_filename = "tmp1234.dataset"
        patch_verify.return_value = (True, tmp_filename)
        
        exp_id = 3
        patch_start.return_value = exp_id
        
        f = user_upload(request)
        self.assertEqual(request.application_url + "/upload/" + str(exp_id), f.location)
        
        patch_check.assert_called_once_with(request, [])
        patch_verify.assert_called_once_with(request)
        patch_start.assert_called_once_with(field_dict, tmp_filename)
        

    @patch('ptmscout.views.experiment.upload.verify_datafile')
    @patch('ptmscout.views.experiment.upload.check_required_fields')
    def test_view_should_handle_form_submission_fail_on_datafile_failure(self, patch_check, patch_verify):
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.user = createMockUser("username", "email", "password", 1)
        
        field_dict = {'name':"something"}
        patch_check.return_value = (True, None, field_dict)
        patch_verify.return_value = (False, strings.failure_reason_experiment_header_no_peptide_column)
        
        result = user_upload(request)
        
        self.assertEqual(strings.failure_reason_experiment_header_no_peptide_column, result['reason'])
        self.assertEqual(field_dict, json.loads(base64.b64decode(result['formfields'])))
        patch_check.assert_called_once_with(request, [])
        patch_verify.assert_called_once_with(request)
        
        
    @patch('ptmscout.views.experiment.upload.check_required_fields')
    def test_view_should_handle_form_submission_fail_on_check(self, patch_check):
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.user = createMockUser("username", "email", "password", 1)
        
        field_dict = {'name':"something"}
        patch_check.return_value = (False, strings.failure_reason_required_fields_cannot_be_empty, field_dict)

        result = user_upload(request)

        self.assertEqual(strings.failure_reason_required_fields_cannot_be_empty, result['reason'])
        self.assertEqual(field_dict, json.loads(base64.b64decode(result['formfields'])))
        patch_check.assert_called_once_with(request, [])
        
    
    def test_view_should_get_users_experiments(self):
        request = DummyRequest()
        
        request.user = createMockUser("username", "email", "password", 1)
        e1 = createMockExperiment(2, 0, 0)
        e2 = createMockExperiment(3, 0, 0)
        p1 = createMockPermission(request.user, e1, 'owner')
        p2 = createMockPermission(request.user, e2, 'view')
        request.user.permissions = [p1,p2]
        
        result = user_upload(request)
        
        expected_experiments = [{'URL':e1.URL, 'parent_id': 0, 'id': 2, 'name': 'Experiment Name2', 'public': 0}]
                
        self.assertEqual(strings.upload_page_title, result['pageTitle'])
        
        del result['user_experiments'][0]['method_calls']
        decoded_data = json.loads(base64.b64decode(result['json_user_data']))
        del decoded_data[0]['method_calls']
        
        self.assertEqual(expected_experiments, result['user_experiments'])
        self.assertEqual(expected_experiments, decoded_data)
        self.assertEqual({}, json.loads(base64.b64decode(result['formfields'])))
        self.assertEqual(None, result['reason'])
        
        
class IntegrationTestUploadView(IntegrationTestCase):
    def test_view_integration(self):
        self.bot.acquire_experiments([26])
        result = self.ptmscoutapp.get("/upload", status=200)
        result.mustcontain(self.bot.user.permissions[0].experiment.name)