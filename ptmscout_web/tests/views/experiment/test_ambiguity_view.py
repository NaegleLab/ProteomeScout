from pyramid.testing import DummyRequest
from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from mock import patch
from tests.views.mocking import createMockExperiment, createMockMeasurement,\
        createMockAmbiguity, createMockProtein, createMockPeptide,\
        createMockPeptideModification, createMockPTM
from ptmscout.views.experiment import ambiguity_view
from ptmscout.config import strings
from ptmscout.utils import forms

class ExperimentAmbiguityViewIntegrationTests(IntegrationTestCase):
    def test_experiment_ambiguity_view(self):
        result = self.ptmscoutapp.get('/experiments/26/ambiguity')

class ExperimentAmbiguityViewTests(UnitTestCase):

    def test_prepopulate_upload_session(self):
        self.fail("not implemented")

    def test_create_ambiguity_schema(self):
        request = DummyRequest()
        request.POST['submitted'] = 'true'

        p = createMockProtein()
        p.name = "mock protein"
        p.species.name = "homo sapiens"

        m1 = createMockMeasurement(4, 1)
        createMockAmbiguity('P01234', m1)
        createMockAmbiguity('P01234-2', m1)
        createMockAmbiguity('P01234-3', m1)
        m1.protein = p
        pep1 = createMockPeptide(4, site_pos=111, site_type='Y')
        createMockPeptideModification(m1, pep1, createMockPTM(name='phosphorylation'))

        m2 = createMockMeasurement(4, 1)
        createMockAmbiguity('P01111', m2)
        createMockAmbiguity('P01111-2', m2)
        createMockAmbiguity('Q0Z456', m2)
        m2.protein = p

        m3 = createMockMeasurement(4, 1)
        createMockAmbiguity('P05500', m3)
        createMockAmbiguity('P05500-2', m3)
        m3.protein = p

        schema, pep_list = ambiguity_view.create_ambiguity_schema([m1,m2,m3], request)

        self.assertEqual({ 'ms%d' % (m1.id), 'ms%d' % (m2.id), 'ms%d' % (m3.id) }, set(schema.field_names.keys()))
        self.assertEqual([ (a,a) for a in ['P01234', 'P01234-2', 'P01234-3'] ], schema.enum_values['ms%d' % (m1.id)])
        self.assertEqual([ (a,a) for a in ['P01111', 'P01111-2', 'Q0Z456'] ], schema.enum_values['ms%d' % (m2.id)])
        self.assertEqual([ (a,a) for a in ['P05500', 'P05500-2'] ], schema.enum_values['ms%d' % (m3.id)])

        self.assertEqual([(0, m1.id, p.name, p.species.name, m1.peptide, [('Y111', 'phosphorylation')]),
                            (1, m2.id, p.name, p.species.name, m2.peptide, []),
                            (2, m3.id, p.name, p.species.name, m3.peptide, [])], pep_list)

    @patch('ptmscout.views.experiment.ambiguity_view.prepopulate_upload_session')
    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.views.experiment.ambiguity_view.create_ambiguity_schema')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_ambiguity_should_eval_on_submission_fail_if_invalid(self, patch_getExperiment, patch_getPeptides, patch_createSchema, patch_validator, patch_createSession):
        request = DummyRequest()
        request.POST['submitted'] = 'true'
        request.matchdict['id'] = '1323'
        request.user = None

        validator = patch_validator.return_value

        session_id = 3000
        schema = forms.FormSchema()
        pep_list = "Some peptide list"
        exp = createMockExperiment(1323)
        patch_getExperiment.return_value = exp
        patch_getPeptides.return_value = "Some peptides"
        patch_createSchema.return_value = schema, pep_list
        patch_createSession.return_value = session_id
        validator.validate.return_value = []

        try:
            result = ambiguity_view.experiment_ambiguity_view(request)
        except HTTPFound, f:
            self.assertEqual("%s/upload/%d/metadata" % (request.application_url, session_id), f.location)
        else:
            self.fail("Expected HTTPFound")

        patch_createSession.assert_called_once_with(exp, request.user, form_schema)
        patch_getExperiment.assert_called_once_with(1323, user=request.user)
        patch_getPeptides.assert_called_once_with(1323, user=request.user)
        patch_createSchema.assert_called_once_with("Some peptides", request)
        validator.validate.assert_called_once_with()

        self.assertEqual([], result['errors'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(pep_list, result['peptides'])
        self.assertEqual(schema, result['formrenderer'].schema)
        self.assertEqual(strings.experiment_ambiguity_page_title, result['pageTitle'])



    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.views.experiment.ambiguity_view.create_ambiguity_schema')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_ambiguity_should_eval_on_submission_fail_if_invalid(self, patch_getExperiment, patch_getPeptides, patch_createSchema, patch_validator):
        request = DummyRequest()
        request.POST['submitted'] = 'true'
        request.matchdict['id'] = '1323'
        request.user = None

        validator = patch_validator.return_value

        schema = forms.FormSchema()
        pep_list = "Some peptide list"
        exp = createMockExperiment(1323)
        patch_getExperiment.return_value = exp
        patch_getPeptides.return_value = "Some peptides"
        patch_createSchema.return_value = schema, pep_list
        validator.validate.return_value = ["Some error"]

        result = ambiguity_view.experiment_ambiguity_view(request)

        patch_getExperiment.assert_called_once_with(1323, user=request.user)
        patch_getPeptides.assert_called_once_with(1323, user=request.user)
        patch_createSchema.assert_called_once_with("Some peptides", request)
        validator.validate.assert_called_once_with()

        self.assertEqual(["Some error"], result['errors'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(pep_list, result['peptides'])
        self.assertEqual(schema, result['formrenderer'].schema)
        self.assertEqual(strings.experiment_ambiguity_page_title, result['pageTitle'])

    @patch('ptmscout.views.experiment.ambiguity_view.create_ambiguity_schema')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_ambiguity(self, patch_getExperiment, patch_getPeptides, patch_createSchema):
        request = DummyRequest()
        request.matchdict['id'] = '1323'
        request.user = None
   
        schema = forms.FormSchema()
        pep_list = "Some peptide list"
        exp = createMockExperiment(1323)
        patch_getExperiment.return_value = exp
        patch_getPeptides.return_value = "Some peptides"
        patch_createSchema.return_value = schema, pep_list

        result = ambiguity_view.experiment_ambiguity_view(request)

        patch_getExperiment.assert_called_once_with(1323, user=request.user)
        patch_getPeptides.assert_called_once_with(1323, user=request.user)
        patch_createSchema.assert_called_once_with("Some peptides", request)

        self.assertEqual([], result['errors'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(pep_list, result['peptides'])
        self.assertEqual(schema, result['formrenderer'].schema)
        self.assertEqual(strings.experiment_ambiguity_page_title, result['pageTitle'])
