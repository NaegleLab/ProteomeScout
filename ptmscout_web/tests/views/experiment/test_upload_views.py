from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
from ptmscout.views.experiment.upload import user_upload
from ptmscout.config import strings
from tests.views.mocking import createMockUser, createMockPermission,\
    createMockExperiment

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
                
        self.assertEqual(strings.upload_page_title, result['pageTitle'])
        self.assertEqual([(e1.id, e1.name)], result['user_experiments'])
        
        
        
class IntegrationTestUploadView(IntegrationTestCase):
    def test_view_integration(self):
        self.bot.acquire_experiments([26])
        result = self.ptmscoutapp.get("/upload", status=200)
        result.mustcontain(self.bot.user.permissions[0].experiment.name)