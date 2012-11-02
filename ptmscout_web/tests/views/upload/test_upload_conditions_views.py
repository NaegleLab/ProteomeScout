from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.database import upload
from ptmscout.config import strings
from ptmscout.utils import forms
from mock import Mock, patch
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockSession, createMockUser,\
    createMockExperiment
from ptmscout.views.upload.upload_conditions import upload_conditions_view,\
    get_form_schema, MAX_VALUES, CONDITION_TYPES, save_form_data
from pyramid.httpexceptions import HTTPFound

class TestUploadConditionsView(UnitTestCase):
    
    def test_save_form_data_should_create_experiment_objects(self):
        request = DummyRequest()
        
        request.POST['0_type'] = "cell"
        request.POST['0_value'] = "v0"
        
        request.POST['1_type'] = "tissue"
        request.POST['1_value'] = "v1"
        
        request.POST['2_type'] = "drug"
        request.POST['2_value'] = "v2"
        
        request.POST['10_type'] = "environment"
        request.POST['10_value'] = "v10"
        
        schema, added_fields = get_form_schema(request)
        exp = createMockExperiment()
        
        save_form_data(exp, schema, added_fields)
        
        exp.addExperimentCondition.assert_any_call("cell", "v0")
        exp.addExperimentCondition.assert_any_call("tissue", "v1")
        exp.addExperimentCondition.assert_any_call("drug", "v2")
        exp.addExperimentCondition.assert_any_call("environment", "v10")
        
        exp.saveExperiment.assert_called_once_with()
        
    
    def test_get_form_schema_should_create_form_and_set_values(self):
        request = DummyRequest()
        
        request.POST['submitted'] = "true"
        request.POST['0_type'] = "f0"
        request.POST['0_value'] = "v0"
        request.POST['1_value'] = "v1"
        request.POST['2_type'] = "f2"
        request.POST['2_value'] = "v2"
        request.POST['10_value'] = "v10"
        
        schema, added_fields = get_form_schema(request)
        
        
        for i in xrange(0, MAX_VALUES):
            self.assertIn('%d_type' % (i), schema.enum_fields)
            self.assertEqual(CONDITION_TYPES, schema.enum_values['%d_type' % (i)])
            self.assertEqual(forms.FormSchema.SELECT, schema.field_types['%d_type' % (i)])
            self.assertEqual(forms.FormSchema.TEXT, schema.field_types['%d_value' % (i)])
            
            parent, _ = schema.conditional_fields['%d_value' % (i)]
            self.assertEqual('%d_type' % (i),parent)
            
        self.assertEqual("f0", schema.form_values['0_type'])
        self.assertEqual("v0", schema.form_values['0_value'])
        self.assertEqual(None, schema.form_values['1_type'])
        self.assertEqual("v1", schema.form_values['1_value'])
        self.assertEqual("f2", schema.form_values['2_type'])
        self.assertEqual("v2", schema.form_values['2_value'])
        self.assertEqual("v10", schema.form_values['10_value'])
        
        self.assertEqual(MAX_VALUES * 2, len(schema.field_names))
        self.assertEqual(set([0,1,2,10]), added_fields)
        
    
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.views.upload.upload_conditions.save_form_data')
    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.database.upload.getSessionById')
    @patch('ptmscout.views.upload.upload_conditions.get_form_schema')
    def test_view_should_validate_on_submission_store_experiment_conditions(self, patch_getSchema, patch_getSession, patch_validator, patch_save_form_data, patch_getExperiment):
        mockValidator = patch_validator.return_value
        mockValidator.validate.return_value = [] 
        
        mockSchema = Mock(spec=forms.FormSchema)
        patch_getSchema.return_value = mockSchema, set([1,2,3])
        
        user = createMockUser()
        exp = createMockExperiment()
        patch_getExperiment.return_value = exp
        
        session = createMockSession(user, experiment_id=exp.id)

        patch_getSession.return_value = session
        request = DummyRequest()
        request.user = user
        request.matchdict['id'] = session.id
        request.POST['submitted'] = "true"
        
        result = upload_conditions_view(request)
        
        self.assertTrue(isinstance(result, HTTPFound))
        self.assertEqual(request.application_url + "/upload/%d/confirm" % (session.id), result.location)
        
        patch_getExperiment.assert_called_once_with(exp.id, request.user, False)
        mockValidator.validate.assert_called_once_with()
        
        patch_save_form_data.assert_called_once_with(exp, mockSchema, set([1,2,3]))
        patch_getSchema.assert_called_once_with(request)
        patch_getSession.assert_called_once_with(session.id, user)
        
    
    @patch('ptmscout.utils.forms.FormValidator')    
    @patch('ptmscout.database.upload.getSessionById')
    @patch('ptmscout.views.upload.upload_conditions.get_form_schema')
    def test_view_should_validate_on_submission_show_errors(self, patch_getSchema, patch_getSession, patch_validator):
        mockValidator = patch_validator.return_value
        mockValidator.validate.return_value = ["some errors"]
        
        mockSchema = Mock(spec=forms.FormSchema)
        patch_getSchema.return_value = mockSchema, set([1,2,3])
        
        user = createMockUser()
        session = createMockSession(user)
        
        patch_getSession.return_value = session
        request = DummyRequest()
        request.user = user
        request.matchdict['id'] = session.id
        request.POST['submitted'] = "true"
        
        result = upload_conditions_view(request)
        
        mockValidator.validate.assert_called_once_with()
        patch_getSchema.assert_called_once_with(request)
        patch_getSession.assert_called_once_with(session.id, user)
        
        self.assertEqual(set([1,2,3]), result['added_fields'])
        self.assertEqual(["some errors"], result['errors'])
        self.assertEqual(strings.experiment_upload_conditions_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_conditions_page_title, result['header'])
        self.assertEqual(mockSchema, result['formrenderer'].schema)
        
        
    
    @patch('ptmscout.database.upload.getSessionById')
    @patch('ptmscout.views.upload.upload_conditions.get_form_schema')    
    def test_view_should_show_form(self, patch_getSchema, patch_getSession):
        mockSchema = Mock(spec=forms.FormSchema)
        patch_getSchema.return_value = mockSchema, set()
        
        user = createMockUser()
        session = createMockSession(user)
        
        patch_getSession.return_value = session
        request = DummyRequest()
        request.user = user
        request.matchdict['id'] = session.id
        
        result = upload_conditions_view(request)
        
        patch_getSchema.assert_called_once_with(request)
        patch_getSession.assert_called_once_with(session.id, user)
        
        self.assertEqual(set(), result['added_fields'])
        self.assertEqual([], result['errors'])
        self.assertEqual(MAX_VALUES, result['MAX_FIELDS'])
        self.assertEqual(strings.experiment_upload_conditions_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_conditions_page_title, result['header'])
        self.assertEqual(mockSchema, result['formrenderer'].schema)
        
        
        
        
class IntegrationTestUploadConditionsView(IntegrationTestCase):
    def test_view_integration(self):
        from ptmscout.database import experiment
        self.bot.login()
        
        exp = experiment.getExperimentById(1, 0, False)
        exp.status = 'preload'
        exp.saveExperiment()
        
        session = upload.Session()
        session.experiment_id = 1
        session.user_id = self.bot.user.id
        session.load_type='new'
        session.data_file='exp_file'
        session.parent_experiment=None
        session.change_description=''
        session.stage='confirm'
        session.save()
        
        result = self.ptmscoutapp.get("/upload/%d/conditions" % session.id, status=200)
        result.mustcontain(strings.experiment_upload_conditions_page_title)
        
        result.form.set('0_type', 'cell')
        result.form.set('0_value', 'HELLA')
        
        result.form.submit(status=302).follow(status=200)
        
        
