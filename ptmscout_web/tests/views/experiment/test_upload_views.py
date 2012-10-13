from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
from ptmscout.views.experiment.upload import user_upload, check_required_fields,\
    verify_datafile, save_datafile, create_experiment_and_start_upload
from ptmscout.config import strings, settings
from tests.views.mocking import createMockUser, createMockPermission,\
    createMockExperiment
import base64
import json
from mock import patch
import os
import cgi

class TestUploadView(UnitTestCase):
    
    def create_field_dict(self):
        field_dict = {}
        field_dict['experiment_name'] = "some name for an exp"
        field_dict['description'] = "Some description"
        field_dict['load_type'] = "extension"
        field_dict['published'] = "yes"
        field_dict['author_contact'] = "someemail@institute.edu"
        field_dict['ambiguous'] = "no"
        
        field_dict['parent_experiment'] = "24"
        field_dict['change_description'] = "We did some stuff to remove some other stuff that wasn't all that great for our analysis" 
        field_dict['authors'] = "Makalaster, David; Claypool, Les; Brigade, Frog"
        field_dict['journal'] = "Journal of serendipitous results"
        field_dict['publication_year'] = "2008"
        field_dict['publication_month'] = "may"
        field_dict['volume'] = "102"
        field_dict['page_start'] = "D102"
        field_dict['page_end'] = ""
        field_dict['pmid'] = "1234"
        field_dict['URL'] = ""
        return field_dict
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_create_experiment_and_start_upload_should_fetch_existing_experiment_for_reload(self, patch_getExperiment):
        field_dict = self.create_field_dict()
        field_dict['load_type'] = 'reload'
        exp_file = "tmpfile1002344.tmp"
        
        current_user = createMockUser("username", "email", "password", "1")
        
        experiment_instance = createMockExperiment(int(field_dict['parent_experiment']), 0, 0)
        patch_getExperiment.return_value = experiment_instance 
        
        result = create_experiment_and_start_upload(field_dict, exp_file, current_user)
        
        patch_getExperiment.assert_called_once_with(int(field_dict['parent_experiment']), current_user)
        experiment_instance.saveExperiment.assert_called_once_with()
        
        self.assertEqual('preload', experiment_instance.status)
        self.assertEqual(experiment_instance.id, result)
        self.assertEqual(field_dict['experiment_name'], experiment_instance.name)
        self.assertEqual(field_dict['author_contact'], experiment_instance.contact)
        self.assertEqual(field_dict['authors'], experiment_instance.author)
        self.assertEqual(field_dict['journal'], experiment_instance.journal)
        self.assertEqual(1, experiment_instance.published)
        self.assertEqual(0, experiment_instance.ambiguity)
        self.assertEqual(0, experiment_instance.public)
        self.assertEqual(1, experiment_instance.export)
        self.assertEqual(int(field_dict['publication_year']), experiment_instance.publication_year)
        self.assertEqual(field_dict['publication_month'], experiment_instance.publication_month)
        self.assertEqual(int(field_dict['volume']), experiment_instance.volume)
        self.assertEqual(field_dict['page_start'], experiment_instance.page_start)
        self.assertEqual(field_dict['page_end'], experiment_instance.page_end)
        self.assertEqual(int(field_dict['pmid']), experiment_instance.PMID)
        self.assertEqual(field_dict['URL'], experiment_instance.URL)
        
        self.assertEqual(exp_file, experiment_instance.dataset)
        
        self.assertEqual(current_user.id, experiment_instance.submitter_id)
    
    @patch('ptmscout.database.experiment.Experiment')
    def test_create_experiment_and_start_upload_should_leave_fields_blank_when_not_needed_and_set_user_permissions(self, patch_experiment):
        field_dict = self.create_field_dict()
        field_dict['load_type'] = 'new'
        field_dict['published'] = 'no'
        exp_file = "tmpfile1002344.tmp"
        
        current_user = createMockUser("username", "email", "password", "1")
        result = create_experiment_and_start_upload(field_dict, exp_file, current_user)
        
        experiment_instance = patch_experiment.return_value
        
        experiment_instance.grantPermission.assert_called_once_with(current_user, 'owner')
        experiment_instance.saveExperiment.assert_called_once_with()
        
        self.assertEqual('preload', experiment_instance.status)
        self.assertEqual(experiment_instance.id, result)
        self.assertEqual(field_dict['experiment_name'], experiment_instance.name)
        self.assertEqual(None, experiment_instance.contact)
        self.assertEqual(None, experiment_instance.author)
        self.assertEqual(None, experiment_instance.journal)
        self.assertEqual(0, experiment_instance.published)
        self.assertEqual(0, experiment_instance.ambiguity)
        self.assertEqual(0, experiment_instance.public)
        self.assertEqual(1, experiment_instance.export)
        self.assertEqual(None, experiment_instance.publication_year)
        self.assertEqual(None, experiment_instance.publication_month)
        self.assertEqual(None, experiment_instance.volume)
        self.assertEqual(None, experiment_instance.page_start)
        self.assertEqual(None, experiment_instance.page_end)
        self.assertEqual(None, experiment_instance.PMID)
        self.assertEqual("", experiment_instance.URL)
        self.assertEqual(None, experiment_instance.experiment_id)
        
        self.assertEqual(exp_file, experiment_instance.dataset)
        self.assertEqual(current_user.id, experiment_instance.submitter_id)
        
        
        
    @patch('ptmscout.database.experiment.Experiment')
    def test_create_experiment_and_start_upload_should_spawn_job(self, patch_experiment):
        field_dict = self.create_field_dict()
        exp_file = "tmpfile1002344.tmp"
        
        current_user = createMockUser("username", "email", "password", "1")
        result = create_experiment_and_start_upload(field_dict, exp_file, current_user)
        
        experiment_instance = patch_experiment.return_value
        experiment_instance.saveExperiment.assert_called_once_with()
        
        self.assertEqual('preload', experiment_instance.status)
        self.assertEqual(experiment_instance.id, result)
        self.assertEqual(field_dict['experiment_name'], experiment_instance.name)
        self.assertEqual(field_dict['author_contact'], experiment_instance.contact)
        self.assertEqual(field_dict['authors'], experiment_instance.author)
        self.assertEqual(field_dict['journal'], experiment_instance.journal)
        self.assertEqual(1, experiment_instance.published)
        self.assertEqual(0, experiment_instance.ambiguity)
        self.assertEqual(0, experiment_instance.public)
        self.assertEqual(1, experiment_instance.export)
        self.assertEqual(int(field_dict['publication_year']), experiment_instance.publication_year)
        self.assertEqual(field_dict['publication_month'], experiment_instance.publication_month)
        self.assertEqual(int(field_dict['volume']), experiment_instance.volume)
        self.assertEqual(field_dict['page_start'], experiment_instance.page_start)
        self.assertEqual(field_dict['page_end'], experiment_instance.page_end)
        self.assertEqual(int(field_dict['pmid']), experiment_instance.PMID)
        self.assertEqual(field_dict['URL'], experiment_instance.URL)
        self.assertEqual(int(field_dict['parent_experiment']), experiment_instance.experiment_id)
        
        self.assertEqual(exp_file, experiment_instance.dataset)
        
        self.assertEqual(current_user.id, experiment_instance.submitter_id)
        
    
    def exec_test_verify(self, filename, tmp_file):
        cwd = settings.ptmscout_path
        
        os.chdir(cwd)
        os.system("mac2unix -n %s %s" % (os.path.join("tests", "behave", "data", filename), os.path.join(settings.experiment_data_file_path, tmp_file)))
        os.system("dos2unix %s" % (os.path.join(settings.experiment_data_file_path, tmp_file)))
        
        error = verify_datafile(tmp_file)
        
        os.remove(os.path.join(settings.experiment_data_file_path, tmp_file))
        
        return error
    
    def test_verify_datafile_should_fail_on_no_acc(self):
        filename = "datasetLoad_noAcc.txt"
        tmp_file = "tmp_test.conv"
        
        error = self.exec_test_verify(filename, tmp_file)
        
        self.assertEqual(strings.failure_reason_experiment_header_no_acc_column, error)
    
    def test_verify_datafile_should_fail_on_multiple_peptide(self):
        filename = "datasetLoad_moreThanOnePep.txt"
        tmp_file = "tmp_test.conv"
        
        error = self.exec_test_verify(filename, tmp_file)
        
        self.assertEqual(strings.failure_reason_experiment_header_multiple_peptide_column, error)
    
    def test_verify_datafile_should_pass_on_good_header(self):
        filename = "datasetLoad_correctDataset.txt"
        tmp_file = "tmp_test.conv"
        
        error = self.exec_test_verify(filename, tmp_file)

        self.assertEqual(None, error)

    @patch('ptmscout.views.experiment.upload.verify_datafile')
    def test_save_datafile_should_save_experiment_datafile_and_pass(self, patch_verify):
        request = DummyRequest()
        request.POST['data_file'] = cgi.FieldStorage()
        
        cwd = settings.ptmscout_path
        os.chdir(cwd)
        
        filename = "datasetLoad_correctDataset.txt"
        request.POST['data_file'].filename = filename
        request.POST['data_file'].file = open(os.sep.join(["tests","behave","data",filename]), 'rb')
        
        patch_verify.return_value = None
        
        status, exp_file = save_datafile(request)
        
        os.chdir(os.sep.join(["tests", "behave", "data"]))
        os.system("mac2unix -n %s tmp_test.conv" % (filename))
        os.system("dos2unix tmp_test.conv")

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
        self.assertEqual(strings.failure_reason_required_fields_cannot_be_empty % "Published", error)
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
    @patch('ptmscout.views.experiment.upload.save_datafile')
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
        patch_start.assert_called_once_with(field_dict, tmp_filename, request.user)
        

    @patch('ptmscout.views.experiment.upload.save_datafile')
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
        
        expected_experiments = [{'status':'loaded', 'URL':e1.URL, 'parent_id': 0, 'id': 2, 'name': 'Experiment Name2', 'public': 0}]
                
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