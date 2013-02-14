from pyramid.testing import DummyRequest
from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from mock import patch, Mock
from tests.views.mocking import createMockExperiment, createMockMeasurement,\
        createMockAmbiguity, createMockProtein, createMockPeptide,\
        createMockPeptideModification, createMockPTM
from ptmscout.views.experiment import comparison_view
from ptmscout.config import strings
from ptmscout.utils import forms
from pyramid.httpexceptions import HTTPFound, HTTPForbidden
from ptmscout.views.experiment import comparison_view


class ComparisonViewIntegrationTests(IntegrationTestCase):
    def test_comparison_view_integration(self):
        result = self.ptmscoutapp.get('/experiments/26/compare')


class ComparisonViewTests(UnitTestCase):
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_comparison_view(self, patch_getExp):
        request = DummyRequest()
        request.matchdict['id'] = '1302'
        request.user=None
        exp = createMockExperiment()
        patch_getExp.return_value = exp

        result = comparison_view.experiment_comparison_view(request)


        self.assertEqual(exp, result['experiment'])
        self.assertEqual([], result['results'])
        self.assertEqual(strings.experiment_compare_page_title, result['pageTitle'])
        patch_getExp.assert_called_once_with(1302, user=request.user)
