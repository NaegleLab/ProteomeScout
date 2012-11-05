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
        
        self.form_schema.add_textarea_field('description', 'Description', 50, 5)
        self.form_schema.add_radio_field('load_type', 'Load Type', [('new', 'New'),('reload', 'Reload')], 'new')
        
        self.form_schema.set_required_field('published')
        self.form_schema.set_field_required_condition('publication_year', 'published', forms.field_not_empty_test)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def test_form_renderer(self):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['publication_year'] = ""
        request.POST['description'] = "stuff"
        request.POST['load_type'] = "reload"
        
        self.form_schema.parse_fields(request)
        renderer = forms.FormRenderer(self.form_schema)
        
        self.assertEqual('<input type="text" id="pmid" class="text_field" size="11" maxlength="10" name="pmid" value="1234" />', renderer.render('pmid', 'pmid', 'text_field').__html__())
        self.assertEqual('<input type="text" id="url" class="text_field" size="30" maxlength="100" name="URL" value="http://somestuff.com" />', renderer.render('URL', 'url', 'text_field').__html__())
        self.assertEqual('<select   name="published">\n<option value="" selected></option>\n<option value="no" >No</option>\n<option value="yes" >Yes</option>\n</select>', renderer.render('published').__html__())
        self.assertEqual('<input type="text"   size="5" maxlength="4" name="publication_year" value="" />', renderer.render('publication_year').__html__())
        self.assertEqual('<textarea cols="50" rows="5" name="description">stuff</textarea>', renderer.render('description').__html__())
        
        new_radio = renderer.render('load_type')
        reload_radio = renderer.render('load_type')
        
        self.assertEqual('<input type="radio"   name="load_type" value="new"  /> New', new_radio.__html__())
        self.assertEqual('<input type="radio"   name="load_type" value="reload" checked /> Reload', reload_radio.__html__())

    def test_check_required_fields_should_fail_on_enum_conflict(self):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "blah"
        request.POST['publication_year'] = "2008"
        
        self.form_schema.parse_fields(request)
        errors = forms.FormValidator(self.form_schema).validate()
        
        self.assertEqual([strings.failure_reason_field_value_not_valid % "Published"], errors)
        self.assertEqual({'pmid':"1234", 'description':None, 'load_type':None, 'URL':"http://somestuff.com", 'published':"blah", 'publication_year':"2008"}, self.form_schema.form_values)
        

    def test_check_required_fields_should_fail_on_numeric_field_format_incorrect(self):
        request = DummyRequest()
        request.POST['pmid'] = "na1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "yes"
        request.POST['publication_year'] = "2008"
        
        self.form_schema.parse_fields(request)
        errors = forms.FormValidator(self.form_schema).validate()
        
        self.assertEqual([strings.failure_reason_field_must_be_numeric % "PubMed ID"], errors)
        self.assertEqual({'pmid':"na1234", 'description':None, 'load_type':None, 'URL':"http://somestuff.com", 'published':"yes", 'publication_year':"2008"}, self.form_schema.form_values)
        

    def test_check_required_fields_should_return_false_if_empty(self):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "   "
        
        self.form_schema.parse_fields(request)
        errors = forms.FormValidator(self.form_schema).validate()
        
        self.assertEqual([strings.failure_reason_required_fields_cannot_be_empty % "Published"], errors)
        self.assertEqual({'pmid':"1234", 'description':None, 'load_type':None, 'publication_year':None, 'URL':"http://somestuff.com", 'published':""}, self.form_schema.form_values)

    def test_check_required_fields_should_return_true_on_success(self):
        request = DummyRequest()
        request.POST['pmid'] = "1234"
        request.POST['URL'] = "http://somestuff.com"
        request.POST['published'] = "yes"
        request.POST['publication_year'] = "2008"
        
        self.form_schema.parse_fields(request)
        errors = forms.FormValidator(self.form_schema).validate()
        
        self.assertEqual([], errors)
        self.assertEqual({'pmid':"1234", 'description':None, 'load_type':None, 'URL':"http://somestuff.com", 'published':"yes", 'publication_year':"2008"}, self.form_schema.form_values)
        
