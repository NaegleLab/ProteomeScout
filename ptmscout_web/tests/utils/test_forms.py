import unittest
from ptmscout.config import strings
from pyramid.testing import DummyRequest
from ptmscout.utils import forms


class FormUtilsTestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
        
        self.form_schema = forms.FormSchema()
        
        self.form_schema.add_numeric_field('pmid', 'PubMed ID', 10)
        self.form_schema.add_text_field('URL', 'URL', 100, 30)
        self.form_schema.add_select_field('published', 'Published', [('no', 'No'),('yes', 'Yes')])
        self.form_schema.add_numeric_field('publication_year', 'Publication Year', 4)
        self.form_schema.add_password_field('password', 'Password', 20)
        
        self.form_schema.add_textarea_field('description', 'Description', 50, 5)
        self.form_schema.add_radio_field('load_type', 'Load Type', [('new', 'New'),('reload', 'Reload')], 'new')

        self.form_schema.add_tooltip('pmid', "this is a tooltip")
        self.form_schema.add_tooltip('description', "this is a textarea tooltip")
        
        self.form_schema.set_required_field('published')
        self.form_schema.set_field_required_condition('publication_year', 'published', forms.field_not_empty_test)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        

    def greater_than_1980(self, field_name, field_value, schema):
        if int(field_value) < 1980:
            return "Error: field '%s' required to be after 1980" % (field_name)
       
    def test_parse_fields_dict_with_non_string_elements_should_convert(self):
        field_dict = {}
        field_dict['pmid'] = 2000
        field_dict['URL'] = "http://somestuff.com"
        field_dict['published'] = "yes"
        field_dict['publication_year'] = 1979
        field_dict['description'] = "stuff"
        field_dict['load_type'] = "reload"

        self.form_schema.parse_fields_dict(field_dict)

        self.assertEqual("2000", self.form_schema.get_form_value('pmid'))
        self.assertEqual("http://somestuff.com", self.form_schema.get_form_value('URL'))
        self.assertEqual("yes", self.form_schema.get_form_value('published'))
        self.assertEqual("1979", self.form_schema.get_form_value('publication_year'))
        self.assertEqual("stuff", self.form_schema.get_form_value('description'))
        self.assertEqual("reload", self.form_schema.get_form_value('load_type'))
 
    def test_add_alternate_validator_should_call_extra_validator(self):
        request = DummyRequest()
        request.POST['pmid'] = "2000"
        request.GET['URL'] = "http://somestuff.com"
        request.POST['published'] = "yes"
        request.POST['publication_year'] = "1979"
        request.GET['description'] = "stuff"
        request.POST['load_type'] = "reload"
        
        self.form_schema.parse_fields(request)
        
        validator = forms.FormValidator(self.form_schema)
        validator.add_alternate_validator('publication_year', self.greater_than_1980)
        validator.add_alternate_validator('pmid', self.greater_than_1980)
        
        errors = validator.validate()
        
        self.assertEqual("Error: field 'Publication Year' required to be after 1980", errors[0])
        
        
        
        
    def test_empty_select_field_required_should_validate_with_error_with_condition_met(self):
        schema = forms.FormSchema()
        
        parent_experiment_options = [(str(24), 'Experiment 1'), (str(28), 'Experiment 2')]
        schema.add_select_field('parent_experiment', 'Parent Experiment', parent_experiment_options)
        schema.add_radio_field('load_type', "Load Type", [('new',"New"),('append',"Append"),('reload',"Reload"),('extension',"Extension")])
        
        schema.set_field_required_condition('parent_experiment', 'load_type', lambda pval: pval != "new")
        
        request = DummyRequest()
        request.POST['load_type'] = 'reload'
        request.POST['parent_experiment'] = ''
        
        schema.parse_fields(request)
        
        errors = forms.FormValidator(schema).validate()
        
        self.assertEqual([strings.failure_reason_required_fields_cannot_be_empty % ('Parent Experiment')], errors)
        
    def test_form_renderer(self):
        request = DummyRequest()
        request.POST = {}
        request.GET = {}
        
        request.POST['pmid'] = "1234"
        request.GET['URL'] = "http://somestuff.com"
        request.POST['publication_year'] = ""
        request.GET['description'] = "stuff"
        request.POST['load_type'] = "reload"
        request.POST['password'] = "somepass"
        
        self.form_schema.parse_fields(request)
        renderer = forms.FormRenderer(self.form_schema)
        
        self.assertEqual('<input type="text" id="pmid" class="text_field" title="this is a tooltip" size="11" maxlength="10" name="pmid" value="1234" />', renderer.render('pmid', 'pmid', 'text_field').__html__())
        self.assertEqual('<input type="text" id="url" class="text_field"  size="30" maxlength="100" name="URL" value="http://somestuff.com" />', renderer.render('URL', 'url', 'text_field').__html__())
        self.assertEqual('<select    name="published">\n<option value="" selected></option>\n<option value="no" >No</option>\n<option value="yes" >Yes</option>\n</select>', renderer.render('published').__html__())
        self.assertEqual('<input type="text"    size="5" maxlength="4" name="publication_year" value="" />', renderer.render('publication_year').__html__())
        self.assertEqual('<input type="password"     maxlength="20" name="password" />', renderer.render('password').__html__())
        self.assertEqual('<textarea id="text-area"  title="this is a textarea tooltip" cols="50" rows="5" name="description">stuff</textarea>', renderer.render('description', 'text-area').__html__())
        
        new_radio = renderer.render('load_type')
        reload_radio = renderer.render('load_type')
        
        self.assertEqual('<input type="radio"    name="load_type" value="new"  /> New', new_radio.__html__())
        self.assertEqual('<input type="radio"    name="load_type" value="reload" checked /> Reload', reload_radio.__html__())

    def test_check_required_fields_should_fail_on_enum_conflict(self):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "blah"
        request.POST['publication_year'] = "2008"
        
        self.form_schema.parse_fields(request)
        errors = forms.FormValidator(self.form_schema).validate()
        
        self.assertEqual([strings.failure_reason_field_value_not_valid % "Published"], errors)
        self.assertEqual({'pmid':"1234", 'password':None, 'description':None, 'load_type':None, 'URL':"http://somestuff.com", 'published':"blah", 'publication_year':"2008"}, self.form_schema.form_values)
        

    def test_check_required_fields_should_fail_on_numeric_field_format_incorrect(self):
        request = DummyRequest()
        request.POST['pmid'] = "na1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "yes"
        request.POST['publication_year'] = "2008"
        
        self.form_schema.parse_fields(request)
        errors = forms.FormValidator(self.form_schema).validate()
        
        self.assertEqual([strings.failure_reason_field_must_be_numeric % "PubMed ID"], errors)
        self.assertEqual({'pmid':"na1234", 'password':None, 'description':None, 'load_type':None, 'URL':"http://somestuff.com", 'published':"yes", 'publication_year':"2008"}, self.form_schema.form_values)
        

    def test_check_required_fields_should_return_false_if_empty(self):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "   "
        
        self.form_schema.parse_fields(request)
        errors = forms.FormValidator(self.form_schema).validate()
        
        self.assertEqual([strings.failure_reason_required_fields_cannot_be_empty % "Published"], errors)
        self.assertEqual({'pmid':"1234", 'password':None, 'description':None, 'load_type':None, 'publication_year':None, 'URL':"http://somestuff.com", 'published':""}, self.form_schema.form_values)

    def test_check_required_fields_should_return_true_on_success(self):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "yes"
        request.POST['publication_year'] = "2008"
        
        self.form_schema.parse_fields(request)
        errors = forms.FormValidator(self.form_schema).validate()
        
        self.assertEqual([], errors)
        self.assertEqual({'pmid':"1234", 'password':None, 'description':None, 'load_type':None, 'URL':"http://somestuff.com", 'published':"yes", 'publication_year':"2008"}, self.form_schema.form_values)
        
