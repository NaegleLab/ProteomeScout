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
    def test_comparison_view_submitted_selection(self, patch_compare):
        request = DummyRequest()
        request.POST = MultiDict()
        request.POST.add('submitted', 'subset')
        request.POST.add('experiment', '1302')
        request.POST.add('experiment', '1304')
        request.POST.add('experiment', '1356')

        request.user=None
        exp = createMockExperiment()
        patch_compare.return_value = {"Expected":"Results"}
        experiment_list = set([1302,1304,1356])

        result = comparison_view.internal_experiment_comparison_view(request, exp)

        patch_compare.assert_called_once_with(exp, request.user, experiment_list)

        self.assertEqual({"Expected":"Results"}, result['results'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.experiment_compare_page_title % (exp.name), result['pageTitle'])

    @patch('ptmscout.views.experiment.comparison_view.compare_to_all')
    def test_comparison_view_submitted_all(self, patch_compare):
        request = DummyRequest()
        request.POST['submitted'] = 'all'
        request.user=None
        
        exp = createMockExperiment()
        patch_compare.return_value = {"Expected":"Results"}

        result = comparison_view.internal_experiment_comparison_view(request, exp)

        patch_compare.assert_called_once_with(exp, request.user)

        self.assertEqual({"Expected":"Results"}, result['results'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.experiment_compare_page_title % (exp.name), result['pageTitle'])

    def test_comparison_view(self):
        request = DummyRequest()
        request.user=None
        
        exp = createMockExperiment()

        result = comparison_view.internal_experiment_comparison_view(request, exp)


        self.assertEqual(exp, result['experiment'])
        self.assertEqual(None, result['results'])
        self.assertEqual(strings.experiment_compare_page_title % (exp.name), result['pageTitle'])
