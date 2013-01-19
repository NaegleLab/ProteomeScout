import unittest
from pyramid.testing import DummyRequest
from ptmscout.utils import forms, paginate


class FormUtilsTestCase(unittest.TestCase):
    def setUp(self):
        self.form_schema = forms.FormSchema()
        self.form_schema.add_text_field('acc_search', 'Accession', 100, 30)
        self.form_schema.add_text_field('pep_search', 'Peptide', 100, 30)
        self.form_schema.add_select_field('species', 'Species', [('homo sapiens', 'Homo sapiens'),('rattus rattus', 'Rattus rattus')])

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def set_form_fields(self):
        field_dict = {}
        field_dict['acc_search'] = "P03456"
        field_dict['pep_search'] = ""
        field_dict['species'] = 'homo sapiens'
        self.form_schema.parse_fields_dict(field_dict)
        return field_dict

    def test_paginator_should_build_wo_data(self):
        request = DummyRequest()

        pager = paginate.Paginator(self.form_schema, 20)
        pager.parse_parameters(request)
        pager.set_result_size(1001)

        pager_html = pager.build()

        for i in xrange(1, 52):
            if i != 1:
                self.assertNotEqual(-1, pager_html.find("%d</a>" % (i)))
            else:
                self.assertNotEqual(-1, pager_html.find("1</span>"))

        self.assertEqual(-1, pager_html.find("&lt;&lt;</a>"))
        self.assertNotEqual(-1, pager_html.find("&gt;&gt;</a>"))

    def test_paginator_should_build_with_data(self):
        request = DummyRequest()
        request.GET = self.set_form_fields()
        request.GET['page'] = '20'

        pager = paginate.Paginator(self.form_schema, 20)
        pager.parse_parameters(request)
        pager.set_result_size(1001)

        pager_html = pager.build()

        for i in xrange(1, 52):
            if i != 20:
                self.assertNotEqual(-1, pager_html.find("%d</a>" % (i)))
            else:
                self.assertNotEqual(-1, pager_html.find("20</span>"))
        self.assertNotEqual(-1, pager_html.find("&lt;&lt;</a>"))
        self.assertNotEqual(-1, pager_html.find("&gt;&gt;</a>"))
