from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase
from pyramid.testing import DummyRequest
from ptmscout.views.upload.upload_metadata import upload_metadata, check_required_fields,\
    create_experiment_and_mark_status, create_experiment
from ptmscout.config import strings
from tests.views.mocking import createMockUser, createMockExperiment,\
    createMockSession
from mock import patch
from ptmscout.database import upload

class TestUploadView(UnitTestCase):
    
    def create_field_dict(self):
        field_dict = {}
        field_dict['experiment_name'] = "some name for an exp"
        field_dict['description'] = "Some description"
        field_dict['published'] = "yes"
        field_dict['author_contact'] = "someemail@institute.edu"
        field_dict['ambiguous'] = "no"
        
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
    def test_create_experiment_should_fetch_existing_experiment_for_reload(self, patch_getExperiment):
        field_dict = self.create_field_dict()
        
        parent_exp_id=20
        current_user = createMockUser()
        
        session = createMockSession(current_user, data_file="tmpfile1002344.tmp", load_type='reload', parent_experiment=parent_exp_id, stage='metadata')
        
        experiment_instance = createMockExperiment(parent_exp_id)
        patch_getExperiment.return_value = experiment_instance 
        
        result = create_experiment(field_dict, session, current_user)
        
        patch_getExperiment.assert_called_once_with(parent_exp_id, current_user)
        experiment_instance.saveExperiment.assert_called_once_with()

        
        self.assertEqual('preload', experiment_instance.status)
        self.assertEqual(experiment_instance, result)
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
        self.assertEqual(session.data_file, experiment_instance.dataset)
        
        self.assertEqual(current_user.id, experiment_instance.submitter_id)
    
    @patch('ptmscout.database.experiment.Experiment')
    def test_create_experiment_should_leave_fields_blank_when_not_needed_and_set_user_permissions(self, patch_experiment):
        field_dict = self.create_field_dict()
        field_dict['published'] = 'no'
        
        current_user = createMockUser()
        session = createMockSession(current_user, data_file="tmpfile1002344.tmp", load_type='new', parent_experiment=20, change_description='some description', stage='metadata')
        result = create_experiment(field_dict, session, current_user)
        
        experiment_instance = patch_experiment.return_value
        
        experiment_instance.grantPermission.assert_called_once_with(current_user, 'owner')
        experiment_instance.saveExperiment.assert_called_once_with()
        
        self.assertEqual('preload', experiment_instance.status)
        self.assertEqual(experiment_instance, result)
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
        
        self.assertEqual(session.data_file, experiment_instance.dataset)
        self.assertEqual(current_user.id, experiment_instance.submitter_id)
        
        
    @patch('ptmscout.database.experiment.Experiment')
    def test_create_experiment(self, patch_experiment):
        field_dict = self.create_field_dict()
        exp_file = "tmpfile1002344.tmp"
        
        current_user = createMockUser()
        session = createMockSession(current_user, data_file="tmpfile1002344.tmp", load_type='extension', parent_experiment=20, change_description='some description', stage='metadata')
        result = create_experiment(field_dict, session, current_user)
        
        experiment_instance = patch_experiment.return_value
        experiment_instance.saveExperiment.assert_called_once_with()
        
        self.assertEqual('preload', experiment_instance.status)
        self.assertEqual(experiment_instance, result)
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
        self.assertEqual(session.parent_experiment, experiment_instance.experiment_id)
        
        self.assertEqual(exp_file, experiment_instance.dataset)
        
        self.assertEqual(current_user.id, experiment_instance.submitter_id)
    
    @patch('ptmscout.views.upload.upload_metadata.create_experiment')
    def test_create_experiment_and_mark_status(self, patch_create_exp):
        user = createMockUser()
        session = createMockSession(user)
        exp = createMockExperiment()
        
        patch_create_exp.return_value = exp
        
        field_dict= {"some":"data"}
        create_experiment_and_mark_status(field_dict, session, user)
        
        patch_create_exp.assert_called_once_with(field_dict, session, user)
        
        self.assertEqual('confirm', session.stage)
        self.assertEqual(exp.id, session.experiment_id)
        
        session.save.assert_called_once_with()
        
        
        
    @patch('ptmscout.views.upload.upload_metadata.get_enum_fields')
    @patch('ptmscout.views.upload.upload_metadata.get_numeric_fields')
    @patch('ptmscout.views.upload.upload_metadata.get_required_fields')
    def test_check_required_fields_should_fail_on_enum_conflict(self, patch_required, patch_numeric, patch_enum):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "blah"
        request.POST['publication_year'] = "2008"
        
        patch_required.return_value = ['published']
        patch_numeric.return_value = ['publication_year', 'pmid']
        patch_enum.return_value = {'published':set(["yes","no"])}
        
        status, error, field_dict = check_required_fields(request)
        
        patch_required.assert_called_once_with(request)
        patch_numeric.assert_called_once_with(request)
        patch_enum.assert_called_once_with(request)
        
        self.assertEqual(False, status)
        self.assertEqual(strings.failure_reason_field_value_not_valid % "Published", error)
        self.assertEqual({'pmid':"1234", 'URL':"http://somestuff.com", 'published':"blah", 'publication_year':"2008"}, field_dict)
        

    @patch('ptmscout.views.upload.upload_metadata.get_numeric_fields')
    @patch('ptmscout.views.upload.upload_metadata.get_required_fields')
    def test_check_required_fields_should_fail_on_numeric_field_format_incorrect(self, patch_required, patch_numeric):
        request = DummyRequest()
        request.POST['pmid'] = "na1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "yes"
        request.POST['publication_year'] = "2008"
        
        patch_required.return_value = ['published']
        patch_numeric.return_value = ['publication_year', 'pmid']
        
        status, error, field_dict = check_required_fields(request)
        
        patch_required.assert_called_once_with(request)
        patch_numeric.assert_called_once_with(request)
        
        self.assertEqual(False, status)
        self.assertEqual(strings.failure_reason_field_must_be_numeric % "PubMed ID", error)
        self.assertEqual({'pmid':"na1234", 'URL':"http://somestuff.com", 'published':"yes", 'publication_year':"2008"}, field_dict)
        

    @patch('ptmscout.views.upload.upload_metadata.get_required_fields')
    def test_check_required_fields_should_return_false_if_empty(self, patch_required):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "   "
        
        patch_required.return_value = ['published']
        
        status, error, field_dict = check_required_fields(request)
        
        patch_required.assert_called_once_with(request)
        
        self.assertEqual(False, status)
        self.assertEqual(strings.failure_reason_required_fields_cannot_be_empty % "Published", error)
        self.assertEqual({'pmid':"1234", 'URL':"http://somestuff.com", 'published':""}, field_dict)

    @patch('ptmscout.views.upload.upload_metadata.get_numeric_fields')
    @patch('ptmscout.views.upload.upload_metadata.get_enum_fields')
    @patch('ptmscout.views.upload.upload_metadata.get_required_fields')
    def test_check_required_fields_should_return_true_on_success(self, patch_required, patch_enums, patch_numeric):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "yes"
        request.POST['publication_year'] = "2008"
        
        patch_required.return_value = ['published']
        patch_numeric.return_value = ['publication_year']
        patch_enums.return_value = {}
        
        status, error, field_dict = check_required_fields(request)
        
        patch_required.assert_called_once_with(request)
        patch_numeric.assert_called_once_with(request)
        patch_enums.assert_called_once_with(request)
        
        self.assertEqual(True, status)
        self.assertEqual(None, error)
        self.assertEqual({'pmid':"1234", 'URL':"http://somestuff.com", 'published':"yes", 'publication_year':"2008"}, field_dict)
        

    @patch('ptmscout.database.upload.getSessionById')
    @patch('ptmscout.views.upload.upload_metadata.create_experiment_and_mark_status')
    @patch('ptmscout.views.upload.upload_metadata.check_required_fields')
    def test_view_should_start_upload_on_successful_verification(self, patch_check, patch_start, patch_getSession):
        request = DummyRequest()
        
        session_id = 10
        request.matchdict['id'] = session_id
        
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        field_dict = {'name':"something"}
        patch_check.return_value = (True, None, field_dict)
        
        session = createMockSession(request.user)
        patch_getSession.return_value = session 
        
        f = upload_metadata(request)
        self.assertEqual(request.application_url + "/upload/%d/conditions" % session_id, f.location)
        
        patch_check.assert_called_once_with(request)
        patch_start.assert_called_once_with(field_dict, session, request.user)
        
    @patch('ptmscout.database.upload.getSessionById')
    @patch('ptmscout.views.upload.upload_metadata.check_required_fields')
    def test_view_should_handle_form_submission_fail_on_check(self, patch_check, patch_getSession):
        request = DummyRequest()
        
        session_id = 10
        request.matchdict['id'] = session_id
        
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        field_dict = {'name':"something"}
        patch_check.return_value = (False, strings.failure_reason_required_fields_cannot_be_empty, field_dict)

        session = createMockSession(request.user)
        patch_getSession.return_value = session 
        
        result = upload_metadata(request)

        self.assertEqual(strings.failure_reason_required_fields_cannot_be_empty, result['reason'])
        self.assertEqual(field_dict, result['formfields'])
        patch_check.assert_called_once_with(request)
        
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_populate_field_if_parent_experiment_avaiable(self, patch_getSession, patch_getExperiment):
        request = DummyRequest()
        session_id = 10
        request.matchdict['id'] = session_id
        request.user = createMockUser()
        
        exp = createMockExperiment()
        patch_getExperiment.return_value = exp

        exp.published = 1
        exp.page_start = '200'
        exp.page_end = '300'
        exp.URL = 'someurl'
        exp.contact = 'someguy@someplace.edu'
        exp.author = 'Some G'
        exp.name = 'Some name'
        exp.publication_year = 2008
        exp.publication_month = 'february'
        exp.volume = '3'
        exp.journal = "some journal"
        exp.description = "some description"
        exp.PMID = 1223456
        exp.ambiguity = 1
        
        field_dict = {
                     'pmid':exp.PMID,
                     'publication_year': exp.publication_year,
                     'publication_month': exp.publication_month, 
                     'published': 'yes' if exp.published == 1 else 'no',
                     'ambiguous': 'yes' if exp.ambiguity == 1 else 'no',
                     'author_contact' : exp.contact,
                     'page_start': exp.page_start,
                     'page_end': exp.page_end,
                     'authors': exp.author,
                     'journal': exp.journal,
                     'volume': exp.volume,
                     'description': exp.description,
                     'experiment_name':"",
                     'URL': exp.URL,
                     'notes':""
        }

        session = createMockSession(request.user)
        session.parent_experiment = exp.id
        patch_getSession.return_value = session 
        
        result = upload_metadata(request)

        patch_getExperiment.assert_called_once_with(exp.id, request.user)
        self.assertEqual(field_dict, result['formfields'])
        
class IntegrationTestUploadMetadataView(IntegrationTestCase):
    def test_view_integration(self):
        
        session = upload.Session()
        session.data_file = 'some_filename'
        session.change_description = ''
        session.experiment_id = None
        session.load_type = 'new'
        session.parent_experiment = None
        session.stage = 'metadata'
        session.user_id = self.bot.user.id
        session.save()
        
        self.ptmscoutapp.get("/upload/%d/metadata" % (session.id), status=200)