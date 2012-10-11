from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
from ptmscout.views.experiment.upload import user_upload
from ptmscout.config import strings
from tests.views.mocking import createMockUser, createMockPermission,\
    createMockExperiment
import base64
import json

class TestUploadView(UnitTestCase):
    def test_view_should_get_users_experiments(self):
        request = DummyRequest()
        
        request.user = createMockUser("username", "email", "password", 1)
        e1 = createMockExperiment(2, 0, 0)
        e2 = createMockExperiment(3, 0, 0)
        p1 = createMockPermission(request.user, e1, 'owner')
        p2 = createMockPermission(request.user, e2, 'view')
        request.user.permissions = [p1,p2]
        
        result = user_upload(request)
        
        expected_experiments = [{'URL':'url', 'parent_id': 0, 'id': 2, 'name': 'Experiment Name2', 'public': 0}]
                
        self.assertEqual(strings.upload_page_title, result['pageTitle'])
        
        del result['user_experiments'][0]['method_calls']
        decoded_data = json.loads(base64.b64decode(result['json_user_data']))
        del decoded_data[0]['method_calls']
        
        self.assertEqual(expected_experiments, result['user_experiments'])
        self.assertEqual(expected_experiments, decoded_data)
        
        
class IntegrationTestUploadView(IntegrationTestCase):
    def test_view_integration(self):
        self.bot.acquire_experiments([26])
        result = self.ptmscoutapp.get("/upload", status=200)
        result.mustcontain(self.bot.user.permissions[0].experiment.name)