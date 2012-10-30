from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
import json
from ptmscout.webservice.query import call_get_pubmed_record_by_id
from pyramid.testing import DummyRequest

class IntegrationTestQueryView(IntegrationTestCase):
    def test_get_pubmed_record_should_return_json_format_data(self):
        result = self.ptmscoutapp.get("/webservice/pubmed/12230038", status=200)
        result.json
        
    def test_get_pubmed_should_be_restricted(self):
        self.bot.logout()
        result = self.ptmscoutapp.get("/webservice/pubmed/12230038", status=200)
        result.mustcontain("Forbidden")