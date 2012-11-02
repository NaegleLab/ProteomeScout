from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
import json
from ptmscout.webservice.query import call_get_pubmed_record_by_id,\
    get_autocomplete_for_field
from pyramid.testing import DummyRequest
from mock import patch
from tests.views.mocking import createMockUser

class TestQuery(UnitTestCase):
    @patch('ptmscout.database.experiment.getValuesForField')
    def test_get_autocomplete_for_field(self, patch_getValues):
        request = DummyRequest()
        request.matchdict['field'] = "fieldname"
        request.user = createMockUser()
        
        patch_getValues.return_value = ["some","values"]
        
        result = get_autocomplete_for_field(request)
        
        patch_getValues.assert_called_once_with("fieldname", request.user)
        
        self.assertEqual({"fieldname":["some","values"]}, result)
        
class IntegrationTestQueryView(IntegrationTestCase):
    def test_get_pubmed_record_should_return_json_format_data(self):
        result = self.ptmscoutapp.get("/webservice/pubmed/12230038", status=200)
        result.json
        
    def test_get_pubmed_should_be_restricted(self):
        self.bot.logout()
        result = self.ptmscoutapp.get("/webservice/pubmed/12230038", status=200)
        result.mustcontain("Forbidden")
        
    def test_get_autocomplete_for_field_is_restricted(self):
        self.bot.logout()
        result = self.ptmscoutapp.get("/webservice/autocomplete/environment", status=200)
        result.mustcontain("Forbidden")
        
    def test_get_autocomplete_for_field(self):
        result = self.ptmscoutapp.get("/webservice/autocomplete/environment", status=200)
        result.json