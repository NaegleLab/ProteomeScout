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

class ComparisonViewIntegrationTests(IntegrationTestCase):
    def test_comparison_view_integration(self):
        result = self.ptmscoutapp.get('/experiments/26/compare')


class ComparisonViewTests(UnitTestCase):
    def test_comparison_view(self):
        pass
