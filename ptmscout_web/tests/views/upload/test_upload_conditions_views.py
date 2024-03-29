from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.database import upload
from ptmscout.config import strings
from ptmscout.utils import forms, wizard
from mock import Mock, patch
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockSession, createMockUser,\
    createMockExperiment, createMockCondition
from ptmscout.views.upload.upload_conditions import upload_conditions,\
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
        
        
        exp = createMockExperiment()
        exp.conditions = []
        schema, added_fields = get_form_schema(exp, None, request)
        
        save_form_data(exp, schema, added_fields)
        
        exp.addExperimentCondition.assert_any_call("cell", "v0")
        exp.addExperimentCondition.assert_any_call("tissue", "v1")
        exp.addExperimentCondition.assert_any_call("drug", "v2")
        exp.addExperimentCondition.assert_any_call("environment", "v10")
        
        exp.saveExperiment.assert_called_once_with()
    
    def test_get_form_schema_should_get_values_from_parent_experiment_if_not_submitted_and_experiment_conditions_empty(self):
        request = DummyRequest()
        
        exp = createMockExperiment()
        parent_exp = createMockExperiment()
        c1 = createMockCondition('f0', 'v0')
        c2 = createMockCondition('f1', 'v1')
        c3 = createMockCondition('f2', 'v2')
        exp.conditions = []
        parent_exp.conditions = [c1,c2,c3]
        
        schema, added_fields = get_form_schema(exp, parent_exp, request)
        
        
        for i in xrange(0, MAX_VALUES):
            self.assertIn('%d_type' % (i), schema.enum_fields)
            self.assertEqual(CONDITION_TYPES, schema.enum_values['%d_type' % (i)])
            self.assertEqual(forms.FormSchema.SELECT, schema.field_types['%d_type' % (i)])
            self.assertEqual(forms.FormSchema.TEXT, schema.field_types['%d_value' % (i)])
            
            parent, _ = schema.conditional_fields['%d_value' % (i)]
            self.assertEqual('%d_type' % (i),parent)
            
        self.assertEqual("f0", schema.form_values['0_type'])
        self.assertEqual("v0", schema.form_values['0_value'])
        self.assertEqual("f1", schema.form_values['1_type'])
        self.assertEqual("v1", schema.form_values['1_value'])
        self.assertEqual("f2", schema.form_values['2_type'])
        self.assertEqual("v2", schema.form_values['2_value'])
        
        self.assertEqual(MAX_VALUES * 2, len(schema.field_names))
        self.assertEqual(set([0,1,2]), added_fields)

    def test_get_form_schema_should_get_values_from_experiment_if_not_submitted(self):
        request = DummyRequest()
        
        exp = createMockExperiment()
        parent_exp = createMockExperiment()
        c1 = createMockCondition('f0', 'v0')
        c2 = createMockCondition('f1', 'v1')
        c3 = createMockCondition('f2', 'v2')
        parent_exp.conditions = []
        exp.conditions = [c1,c2,c3]

        schema, added_fields = get_form_schema(exp, parent_exp, request)
        
        
        for i in xrange(0, MAX_VALUES):
            self.assertIn('%d_type' % (i), schema.enum_fields)
            self.assertEqual(CONDITION_TYPES, schema.enum_values['%d_type' % (i)])
            self.assertEqual(forms.FormSchema.SELECT, schema.field_types['%d_type' % (i)])
            self.assertEqual(forms.FormSchema.TEXT, schema.field_types['%d_value' % (i)])
            
            parent, _ = schema.conditional_fields['%d_value' % (i)]
            self.assertEqual('%d_type' % (i),parent)
            
        self.assertEqual("f0", schema.form_values['0_type'])
        self.assertEqual("v0", schema.form_values['0_value'])
        self.assertEqual("f1", schema.form_values['1_type'])
        self.assertEqual("v1", schema.form_values['1_value'])
        self.assertEqual("f2", schema.form_values['2_type'])
        self.assertEqual("v2", schema.form_values['2_value'])
        
        self.assertEqual(MAX_VALUES * 2, len(schema.field_names))
        self.assertEqual(set([0,1,2]), added_fields)
    
    def test_get_form_schema_should_create_form_and_set_values(self):
        request = DummyRequest()
        
        request.POST['submitted'] = "true"
        request.POST['0_type'] = "f0"
        request.POST['0_value'] = "v0"
        request.POST['1_value'] = "v1"
        request.POST['2_type'] = "f2"
        request.POST['2_value'] = "v2"
        request.POST['10_value'] = "v10"
        
        exp = createMockExperiment()
        schema, added_fields = get_form_schema(exp, None, request)
        
        
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
    @patch('ptmscout.views.upload.upload_conditions.get_form_schema')
    @patch('ptmscout.views.upload.upload_conditions.create_nav_wizard')
    def test_view_should_validate_on_submission_store_experiment_conditions(self, patch_navWizard, patch_getSchema, patch_validator, patch_save_form_data, patch_getExperiment):
        mockValidator = patch_validator.return_value
        mockValidator.validate.return_value = [] 
        
        mockSchema = Mock(spec=forms.FormSchema)
        patch_getSchema.return_value = mockSchema, set([1,2,3])
        patch_navWizard.return_value = Mock(spec=wizard.WizardNavigation)
        
        user = createMockUser()
        exp = createMockExperiment()
        patch_getExperiment.return_value = exp
        
        session = createMockSession(user, experiment_id=exp.id)
        session.stage='condition'
        session.parent_experiment = None
        
        request = DummyRequest()
        request.user = user
        request.matchdict['id'] = session.id
        request.POST['submitted'] = "true"
        
        result = upload_conditions(request, session)
        
        self.assertTrue(isinstance(result, HTTPFound))
        self.assertEqual(request.application_url + "/upload/%d/confirm" % (session.id), result.location)
        self.assertEqual('confirm', session.stage)
        
        session.save.assert_called_once_with()
        
        patch_getExperiment.assert_called_once_with(exp.id, request.user, False)
        mockValidator.validate.assert_called_once_with()
        
        patch_save_form_data.assert_called_once_with(exp, mockSchema, set([1,2,3]))
        patch_getSchema.assert_called_once_with(exp, None, request)
        
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.utils.forms.FormValidator')    
    @patch('ptmscout.views.upload.upload_conditions.get_form_schema')
    @patch('ptmscout.views.upload.upload_conditions.create_nav_wizard')
    def test_view_should_validate_on_submission_show_errors(self, patch_navWizard, patch_getSchema, patch_validator, patch_getExperiment):
        mockValidator = patch_validator.return_value
        mockValidator.validate.return_value = ["some errors"]
        
        mockSchema = Mock(spec=forms.FormSchema)
        patch_getSchema.return_value = mockSchema, set([1,2,3])
        patch_navWizard.return_value = Mock(spec=wizard.WizardNavigation)
        
        user = createMockUser()
        exp = createMockExperiment()
        patch_getExperiment.return_value = exp
        session = createMockSession(user, experiment_id=exp.id)
        session.parent_experiment = None
        
        request = DummyRequest()
        request.user = user
        request.matchdict['id'] = session.id
        request.POST['submitted'] = "true"
        
        result = upload_conditions(request, session)
        
        mockValidator.validate.assert_called_once_with()
        patch_getSchema.assert_called_once_with(exp, None, request)
        
        self.assertEqual(set([1,2,3]), result['added_fields'])
        self.assertEqual(["some errors"], result['errors'])
        self.assertEqual(strings.experiment_upload_conditions_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_conditions_page_title, result['header'])
        self.assertEqual(mockSchema, result['formrenderer'].schema)
        self.assertEqual(patch_navWizard.return_value, result['navigation'])
        
        
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.views.upload.upload_conditions.get_form_schema')    
    @patch('ptmscout.views.upload.upload_conditions.create_nav_wizard')
    def test_view_should_show_form(self, patch_navWizard, patch_getSchema, patch_getExperiment):
        mockSchema = Mock(spec=forms.FormSchema)
        patch_getSchema.return_value = mockSchema, set()
        patch_navWizard.return_value = Mock(spec=wizard.WizardNavigation)
        
        user = createMockUser()
        exp = createMockExperiment()
        patch_getExperiment.return_value = exp
        session = createMockSession(user, experiment_id=exp.id)
        session.parent_experiment = None
        
        request = DummyRequest()
        request.user = user
        request.matchdict['id'] = session.id
        
        result = upload_conditions(request, session)
        
        patch_getSchema.assert_called_once_with(exp, None, request)
        
        self.assertEqual(session.id, result['session_id'])
        self.assertEqual(set(), result['added_fields'])
        self.assertEqual([], result['errors'])
        self.assertEqual(MAX_VALUES, result['MAX_FIELDS'])
        self.assertEqual(strings.experiment_upload_conditions_page_title, result['pageTitle'])
        self.assertEqual(strings.experiment_upload_conditions_page_title, result['header'])
        self.assertEqual(mockSchema, result['formrenderer'].schema)
        self.assertEqual(patch_navWizard.return_value, result['navigation'])
        
        
        
        
class IntegrationTestUploadConditionsView(IntegrationTestCase):
    def test_view_integration(self):
        from ptmscout.database import experiment
        self.bot.login()
        
        session = upload.Session()
        session.experiment_id = 1
        session.user_id = self.bot.user.id
        session.resource_type='experiment'
        session.load_type='new'
        session.data_file='exp_file'
        session.parent_experiment=None
        session.change_name=''
        session.change_description=''
        session.stage='conditions'
        session.save()
        
        result = self.ptmscoutapp.get("/upload/%d/conditions" % session.id, status=200)
        result.mustcontain(strings.experiment_upload_conditions_page_title)

        result.form.set('0_type', 'cell')
        result.form.set('0_value', 'HELLA')
        result.form.set('1_type', 'drug')
        result.form.set('1_value', 'doxycycline')
        result.form.set('2_type', 'drug')
        result.form.set('2_value', 'natalizumab')
        
        result.form.submit().follow(status=200)
        
    def test_view_integration_should_load_conditions_from_experiment_parent(self):
        from ptmscout.database import experiment
        self.bot.login()

        session = upload.Session()
        session.experiment_id = 26
        session.user_id = self.bot.user.id
        session.resource_type='experiment'
        session.load_type='new'
        session.data_file='exp_file'
        session.parent_experiment=28
        session.change_name=''
        session.change_description=''
        session.stage='conditions'
        session.save()

        result = self.ptmscoutapp.get("/upload/%d/conditions" % session.id, status=200)
        result.mustcontain(strings.experiment_upload_conditions_page_title)

        self.assertEqual('drug', result.form.get('0_type').value)
        self.assertEqual('dasatinib', result.form.get('0_value').value)


