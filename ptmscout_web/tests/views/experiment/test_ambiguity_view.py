from pyramid.testing import DummyRequest
from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from mock import patch, Mock
from tests.views.mocking import createMockExperiment, createMockMeasurement,\
        createMockAmbiguity, createMockProtein, createMockPeptide,\
        createMockPeptideModification, createMockPTM, createMockUser
from ptmscout.views.experiment import ambiguity_view
from ptmscout.config import strings
from ptmscout.utils import forms
from pyramid.httpexceptions import HTTPForbidden

class ExperimentAmbiguityViewIntegrationTests(IntegrationTestCase):
    def test_experiment_ambiguity_view(self):
        result = self.ptmscoutapp.get('/experiments/26/ambiguity')
        result = result.form.submit()
        result.mustcontain(strings.experiment_ambiguity_error_no_change)
        result = self.ptmscoutapp.get('/experiments/26/ambiguity?defaults=true')
        result = result.form.submit()
        result = result.follow()
        result = result.follow()
        result.mustcontain("Experiment Information")
        result.mustcontain("[Default Assignments] Time-resolved mass spectrometry of tyrosine phosphorylation sites in the epidermal growth factor receptor signaling network reveals dynamic modules.")

class ExperimentAmbiguityViewTests(UnitTestCase):

    @patch('ptmscout.database.upload.Session')
    @patch('ptmscout.database.experiment.Experiment')
    def test_create_session(self, patch_exp, patch_session):
        exp = createMockExperiment()
        exp.public = 1

        new_session = patch_session.return_value
        new_exp = patch_exp.return_value

        new_session.id = 100
        new_exp.id = 1300

        exp_file = 'filename.tmp'
        name_prefix = '[default]'
        columns = [{'type':'accession','label':''},{'type':'peptide','label':''}]
        units = 'time(min)'
        user = createMockUser()

        session_id = ambiguity_view.create_session(exp, user, exp_file, name_prefix, columns, units)

        new_exp.copyData.assert_called_once_with(exp)
        self.assertEqual(0, new_exp.public)
        self.assertEqual(name_prefix + exp.name, new_exp.name)
        new_exp.grantPermission.assert_called_once_with(user, 'owner')
        new_exp.saveExperiment.assert_called_once_with()

        self.assertEqual(new_session.id, session_id)
        self.assertEqual(exp_file, new_session.data_file)
        self.assertEqual('experiment', new_session.resource_type)
        self.assertEqual('extension', new_session.load_type)
        self.assertEqual('metadata', new_session.stage)
        self.assertEqual('accession', new_session.columns[0].type)
        self.assertEqual('', new_session.columns[0].label)
        self.assertEqual('peptide', new_session.columns[1].type)
        self.assertEqual('', new_session.columns[1].label)

        new_session.save.assert_called_once_with()

    def test_create_ambiguity_schema(self):
        request = DummyRequest()
        request.POST['submitted'] = 'true'

        p = createMockProtein()
        p.name = "mock protein"
        p.acc_gene = "AGENE"
        p.species.name = "homo sapiens"

        m1 = createMockMeasurement(4, 1)
        m1.query_accession = 'gi|12345677'
        createMockAmbiguity('AGENE', 'P01234', m1)
        createMockAmbiguity('AGENE', 'P01234-2', m1)
        createMockAmbiguity('AGENE1', 'P01234-3', m1)
        m1.protein = p
        pep1 = createMockPeptide(4, site_pos=111, site_type='Y')
        createMockPeptideModification(m1, pep1, createMockPTM(name='phosphorylation'))

        m2 = createMockMeasurement(4, 1)
        m2.query_accession = 'Q0Z456'
        createMockAmbiguity('AGENE1', 'P01111', m2)
        createMockAmbiguity('AGENE1', 'P01111-2', m2)
        createMockAmbiguity('AGENE1', 'Q0Z456', m2)
        m2.protein = p

        m3 = createMockMeasurement(4, 1)
        m3.query_accession = 'P05500-2'
        createMockAmbiguity('AGENE', 'P05500', m3)
        createMockAmbiguity('AGENE', 'P05500-2', m3)
        m3.protein = p

        schema, pep_list = ambiguity_view.create_ambiguity_schema([m1,m2,m3], request)

        self.assertEqual({ 'ms%d' % (m1.id), 'ms%d' % (m2.id), 'ms%d' % (m3.id) }, set(schema.field_names.keys()))


        exp_values = \
                [ ('gi|12345677', 'AGENE : gi|12345677'), ('P01234', 'AGENE : P01234'),
                  ('P01234-2', 'AGENE : P01234-2'), ('P01234-3', 'AGENE1 : P01234-3') ]
        self.assertEqual('gi|12345677', schema.field_defaults['ms%d' % (m1.id)])
        self.assertEqual(exp_values, schema.enum_values['ms%d' % (m1.id)])

        exp_values = \
                [('P01111', 'AGENE1 : P01111'), ('P01111-2', 'AGENE1 : P01111-2'), ('Q0Z456', 'AGENE1 : Q0Z456')]
        self.assertEqual('Q0Z456', schema.field_defaults['ms%d' % (m2.id)])
        self.assertEqual(exp_values, schema.enum_values['ms%d' % (m2.id)])


        exp_values = \
                [('P05500', 'AGENE : P05500'), ('P05500-2', 'AGENE : P05500-2')]
        self.assertEqual('P05500-2', schema.field_defaults['ms%d' % (m3.id)])
        self.assertEqual(exp_values, schema.enum_values['ms%d' % (m3.id)])

        self.assertEqual([(0, m1.id, m1.query_accession, p.acc_gene, p.name, p.species.name, m1.peptide, [('Y111', 'phosphorylation')]),
                            (1, m2.id, m2.query_accession, p.acc_gene, p.name, p.species.name, m2.peptide, []),
                            (2, m3.id, m3.query_accession, p.acc_gene, p.name, p.species.name, m3.peptide, [])], pep_list)

    @patch('ptmscout.views.experiment.ambiguity_view.prepopulate_upload_session')
    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.views.experiment.ambiguity_view.create_ambiguity_schema')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')
    def test_experiment_ambiguity_should_eval_on_submission_pass_if_valid(self, patch_getPeptides, patch_createSchema, patch_validator, patch_createSession):
        request = DummyRequest()
        request.POST['submitted'] = 'true'
        request.matchdict['id'] = '1323'
        request.user = None

        validator = patch_validator.return_value

        session_id = 3000
        schema = Mock(spec=forms.FormSchema)
        schema.is_defaulted.return_value=False

        pep_list = "Some peptide list"
        exp = createMockExperiment(1323)
        patch_getPeptides.return_value = "Some peptides"
        patch_createSchema.return_value = schema, pep_list
        patch_createSession.return_value = session_id
        validator.validate.return_value = []

        f = ambiguity_view.internal_experiment_ambiguity_view(request, exp)
        self.assertEqual("%s/upload/%d" % (request.application_url, session_id), f.location)

        patch_createSession.assert_called_once_with(exp, request.user, schema, False)
        patch_getPeptides.assert_called_once_with(1323, user=request.user)
        patch_createSchema.assert_called_once_with("Some peptides", request)
        validator.validate.assert_called_once_with()



    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.views.experiment.ambiguity_view.create_ambiguity_schema')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')
    def test_experiment_ambiguity_should_eval_on_submission_fail_if_invalid(self, patch_getPeptides, patch_createSchema, patch_validator):
        request = DummyRequest()
        request.POST['submitted'] = 'true'
        request.matchdict['id'] = '1323'
        request.user = None

        validator = patch_validator.return_value

        schema = forms.FormSchema()
        pep_list = "Some peptide list"
        exp = createMockExperiment(1323)
        patch_getPeptides.return_value = "Some peptides"
        patch_createSchema.return_value = schema, pep_list
        validator.validate.return_value = ["Some error"]

        result = ambiguity_view.internal_experiment_ambiguity_view(request, exp)

        patch_getPeptides.assert_called_once_with(1323, user=request.user)
        patch_createSchema.assert_called_once_with("Some peptides", request)
        validator.validate.assert_called_once_with()

        self.assertEqual(["Some error", strings.experiment_ambiguity_error_no_change], result['errors'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(pep_list, result['peptides'])
        self.assertEqual(schema, result['formrenderer'].schema)
        self.assertEqual("%s: %s" % (strings.experiment_ambiguity_page_title, exp.name), result['pageTitle'])


    @patch('ptmscout.views.experiment.ambiguity_view.assign_defaults')
    @patch('ptmscout.views.experiment.ambiguity_view.create_ambiguity_schema')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')
    def test_experiment_ambiguity_should_assign_defaults_if_get_defaults_is_true(self, patch_getPeptides, patch_createSchema, patch_assign):
        request = DummyRequest()
        request.matchdict['id'] = '1323'
        request.GET['defaults'] = 'true'
        request.user = None
   
        schema = forms.FormSchema()
        pep_list = "Some peptide list"
        exp = createMockExperiment(1323)
        patch_getPeptides.return_value = "Some peptides"
        patch_createSchema.return_value = schema, pep_list
        patch_assign.return_value = ["Some list of ms ids that changed"]

        result = ambiguity_view.internal_experiment_ambiguity_view(request, exp)

        patch_getPeptides.assert_called_once_with(1323, user=request.user)
        patch_createSchema.assert_called_once_with("Some peptides", request)
        patch_assign.assert_called_once_with("Some peptides", schema, request.user)

        self.assertEqual(True, result['assigned_defaults'])
        self.assertEqual(["Some list of ms ids that changed"], result['changed_default'])
        self.assertEqual([], result['errors'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(pep_list, result['peptides'])
        self.assertEqual(schema, result['formrenderer'].schema)
        self.assertEqual("%s: %s" % (strings.experiment_ambiguity_page_title, exp.name), result['pageTitle'])

    def test_experiment_ambiguity_should_throw_forbidden_if_no_ambiguous_assignments(self):
        request = DummyRequest()
        request.matchdict['id'] = '1323'
        request.user = None
   
        schema = forms.FormSchema()
        pep_list = "Some peptide list"
        exp = createMockExperiment(1323)
        exp.ambiguity = 0

        try:
            ambiguity_view.internal_experiment_ambiguity_view(request, exp)
        except HTTPForbidden:
            pass
        else:
            self.fail("Expected HTTPForbidden")

    @patch('ptmscout.views.experiment.ambiguity_view.create_ambiguity_schema')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')
    def test_experiment_ambiguity(self, patch_getPeptides, patch_createSchema):
        request = DummyRequest()
        request.matchdict['id'] = '1323'
        request.user = None
   
        schema = forms.FormSchema()
        pep_list = "Some peptide list"
        exp = createMockExperiment(1323)
        patch_getPeptides.return_value = "Some peptides"
        patch_createSchema.return_value = schema, pep_list

        result = ambiguity_view.internal_experiment_ambiguity_view(request, exp)

        patch_getPeptides.assert_called_once_with(1323, user=request.user)
        patch_createSchema.assert_called_once_with("Some peptides", request)

        self.assertEqual([], result['errors'])
        self.assertEqual(set(), result['changed_default'])
        self.assertEqual(False, result['assigned_defaults'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(pep_list, result['peptides'])
        self.assertEqual(schema, result['formrenderer'].schema)
        self.assertEqual("%s: %s" % (strings.experiment_ambiguity_page_title, exp.name), result['pageTitle'])
