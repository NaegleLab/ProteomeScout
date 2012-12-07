from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase
from pyramid.testing import DummyRequest
from ptmscout.views.upload.upload_metadata import upload_metadata,\
    write_experiment_properties, create_schema, get_experiment_ref, mark_status
from ptmscout.config import strings
from tests.views.mocking import createMockUser, createMockExperiment,\
    createMockSession
from mock import patch, Mock
from ptmscout.database import upload
from ptmscout.utils import forms

class TestUploadView(UnitTestCase):
    
    @patch('ptmscout.database.experiment.Experiment')
    def test_get_experiment_ref_should_return_new_experiment(self, patch_experiment):
        user = createMockUser()
        
        session = createMockSession(user)
        
        nexp = get_experiment_ref(session, user)
        self.assertEqual(patch_experiment.return_value, nexp)
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_get_experiment_ref_should_return_existing_experiment(self, patch_getExperiment):
        user = createMockUser()
        
        exp = createMockExperiment()
        patch_getExperiment.return_value = exp
        
        session = createMockSession(user)
        session.experiment_id = exp.id
        
        nexp = get_experiment_ref(session, user)
        
        self.assertEqual(exp, nexp)
        
        patch_getExperiment.assert_called_once_with(exp.id, user, False)
    
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
    
    def test_write_experiment_properties_should_leave_fields_blank_when_not_needed(self):
        field_dict = self.create_field_dict()
        field_dict['published'] = 'no'
        
        request = DummyRequest()
        request.POST = field_dict
        
        schema = create_schema(request)
        
        current_user = createMockUser()
        session = createMockSession(current_user, data_file="tmpfile1002344.tmp", load_type='new', parent_experiment=20, change_description='some description', stage='metadata')
        experiment_instance = createMockExperiment()
        write_experiment_properties(experiment_instance, session, schema, current_user)
        
        experiment_instance.grantPermission.assert_called_once_with(current_user, 'owner')
        experiment_instance.saveExperiment.assert_called_once_with()
        
        self.assertEqual('configuration', experiment_instance.status)
        self.assertEqual(field_dict['experiment_name'], experiment_instance.name)
        self.assertEqual(None, experiment_instance.contact)
        self.assertEqual(None, experiment_instance.author)
        self.assertEqual(None, experiment_instance.journal)
        self.assertEqual(0, experiment_instance.published)
        self.assertEqual(0, experiment_instance.ambiguity)
        self.assertEqual(0, experiment_instance.public)
        self.assertEqual(0, experiment_instance.export)
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
        
        
    def test_write_experiment_properties(self):
        field_dict = self.create_field_dict()
        exp_file = "tmpfile1002344.tmp"
        
        request = DummyRequest()
        request.POST = field_dict
        
        schema = create_schema(request)
        
        current_user = createMockUser()
        session = createMockSession(current_user, data_file="tmpfile1002344.tmp", load_type='extension', parent_experiment=20, change_description='some description', stage='metadata')
        experiment_instance = createMockExperiment()
        write_experiment_properties(experiment_instance, session, schema, current_user)
        
        experiment_instance.saveExperiment.assert_called_once_with()
        
        self.assertEqual('configuration', experiment_instance.status)
        self.assertEqual(field_dict['experiment_name'], experiment_instance.name)
        self.assertEqual(field_dict['author_contact'], experiment_instance.contact)
        self.assertEqual(field_dict['authors'], experiment_instance.author)
        self.assertEqual(field_dict['journal'], experiment_instance.journal)
        self.assertEqual(1, experiment_instance.published)
        self.assertEqual(0, experiment_instance.ambiguity)
        self.assertEqual(0, experiment_instance.public)
        self.assertEqual(0, experiment_instance.export)
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
    
    def test_mark_status(self):
        session = createMockSession(createMockUser())
        exp = createMockExperiment(20)
        mark_status(exp, session)
        
        self.assertEqual(20, session.experiment_id)
        self.assertEqual('condition', session.stage)
        session.save.assert_called_once_with()
        

    @patch('ptmscout.views.upload.upload_metadata.mark_status')
    @patch('ptmscout.views.upload.upload_metadata.write_experiment_properties')
    @patch('ptmscout.views.upload.upload_metadata.get_experiment_ref')
    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.views.upload.upload_metadata.create_schema')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_start_upload_on_successful_verification(self, patch_getSession, patch_create_schema, patch_validator, patch_getref, patch_write, patch_mark):
        request = DummyRequest()
        
        schema = Mock(spec=forms.FormSchema)
        patch_create_schema.return_value = schema
        validator = patch_validator.return_value
        
        session_id = 10
        request.matchdict['id'] = session_id
        
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        exp = createMockExperiment()
        patch_getref.return_value = exp
        
        session = createMockSession(request.user)
        patch_getSession.return_value = session 
        
        validator.validate.return_value = []
        
        f = upload_metadata(request)
        self.assertEqual(request.application_url + "/upload/%d/conditions" % session_id, f.location)
                
        patch_getref.assert_called_once_with(session, request.user)
        patch_write.assert_called_once_with(exp, session, schema, request.user)
        patch_mark.assert_called_once_with(exp, session)
        validator.validate.assert_called_once_with()
        
    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.views.upload.upload_metadata.create_schema')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_handle_form_submission_fail_on_check(self, patch_getSession, patch_create_schema, patch_validator):
        request = DummyRequest()
        
        schema = Mock(spec=forms.FormSchema)
        patch_create_schema.return_value = schema
        validator = patch_validator.return_value
        
        session_id = 10
        request.matchdict['id'] = session_id
        
        request.POST['submitted'] = "true"
        request.user = createMockUser()
        
        session = createMockSession(request.user)
        patch_getSession.return_value = session 
        
        validator.validate.return_value = [strings.failure_reason_required_fields_cannot_be_empty]
        
        result = upload_metadata(request)

        self.assertEqual([strings.failure_reason_required_fields_cannot_be_empty], result['errors'])
        self.assertEqual(schema, result['formrenderer'].schema)
        
        validator.validate.assert_called_once_with()
    
    @patch('ptmscout.views.upload.upload_metadata.create_schema')
    @patch('ptmscout.views.upload.upload_metadata.populate_schema_from_experiment')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_populate_field_if_parent_experiment_avaiable(self, patch_getSession, patch_prepopulate, patch_create_schema):
        request = DummyRequest()
        session_id = 10
        request.matchdict['id'] = session_id
        request.user = createMockUser()
        
        schema = Mock(spec=forms.FormSchema)
        patch_create_schema.return_value = schema
        
        session = createMockSession(request.user)
        session.parent_experiment = 100
        patch_getSession.return_value = session 
        
        result = upload_metadata(request)
        
        patch_prepopulate.assert_called_once_with(schema, session, request.user)
        
        self.assertEqual(strings.upload_page_title, result['pageTitle'])
        self.assertEqual([], result['errors'])
        self.assertEqual(schema, result['formrenderer'].schema)

        
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