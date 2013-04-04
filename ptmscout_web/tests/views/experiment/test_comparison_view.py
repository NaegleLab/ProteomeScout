from pyramid.testing import DummyRequest
from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from mock import patch
from tests.views.mocking import createMockExperiment
from ptmscout.views.experiment import comparison_view
from ptmscout.config import strings
from webob.multidict import MultiDict


class ComparisonViewIntegrationTests(IntegrationTestCase):
    def test_comparison_view_integration(self):
        result = self.ptmscoutapp.post('/experiments/26/compare', {'submitted':'all'})
        result.mustcontain("Novel Sites (5)")
        result.mustcontain("Ambiguous Sites (60)")


class ComparisonViewTests(UnitTestCase):

    @patch('ptmscout.views.experiment.comparison_view.compare_to_all')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_comparison_view_submitted_selection(self, patch_getExp, patch_compare):
        request = DummyRequest()
        request.POST = MultiDict()
        request.POST.add('submitted', 'subset')
        request.POST.add('experiment', '1302')
        request.POST.add('experiment', '1304')
        request.POST.add('experiment', '1356')

        request.matchdict['id'] = '1302'

        request.user=None
        exp = createMockExperiment()
        patch_getExp.return_value = exp
        patch_compare.return_value = {"Expected":"Results"}
        experiment_list = set([1302,1304,1356])

        result = comparison_view.experiment_comparison_view(request)

        patch_compare.assert_called_once_with(exp, request.user, experiment_list)

        self.assertEqual({"Expected":"Results"}, result['results'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.experiment_compare_page_title, result['pageTitle'])
        patch_getExp.assert_called_once_with(1302, user=request.user)

    @patch('ptmscout.views.experiment.comparison_view.compare_to_all')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_comparison_view_submitted_all(self, patch_getExp, patch_compare):
        request = DummyRequest()
        request.POST['submitted'] = 'all'
        request.matchdict['id'] = '1302'
        request.user=None
        exp = createMockExperiment()
        patch_getExp.return_value = exp
        patch_compare.return_value = {"Expected":"Results"}

        result = comparison_view.experiment_comparison_view(request)

        patch_compare.assert_called_once_with(exp, request.user)

        self.assertEqual({"Expected":"Results"}, result['results'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.experiment_compare_page_title, result['pageTitle'])
        patch_getExp.assert_called_once_with(1302, user=request.user)

    @patch('ptmscout.database.experiment.getExperimentById')
    def test_comparison_view(self, patch_getExp):
        request = DummyRequest()
        request.matchdict['id'] = '1302'
        request.user=None
        exp = createMockExperiment()
        patch_getExp.return_value = exp

        result = comparison_view.experiment_comparison_view(request)


        self.assertEqual(exp, result['experiment'])
        self.assertEqual(None, result['results'])
        self.assertEqual(strings.experiment_compare_page_title, result['pageTitle'])
        patch_getExp.assert_called_once_with(1302, user=request.user)
