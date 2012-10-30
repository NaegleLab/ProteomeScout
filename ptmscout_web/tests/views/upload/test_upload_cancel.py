from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.views.upload.upload_cancel import cancel_upload_view
from tests.views.mocking import createMockSession, createMockUser
from pyramid.testing import DummyRequest
from ptmscout.config import strings
from mock import patch

class TestCancelView(UnitTestCase):

    @patch('ptmscout.database.upload.getSessionById')
    def test_view(self, patch_session):
        user = createMockUser()
        session = createMockSession(user)
        patch_session.return_value = session  
        
        request = DummyRequest()
        request.matchdict['id'] = '247'
        request.user = user
        
        result = cancel_upload_view(request)

        patch_session.assert_called_once_with(247, user)
        session.delete.assert_called_once_with()
        
        self.assertEqual(strings.cancel_upload_successful_page_title, result['pageTitle'])
        self.assertEqual(strings.cancel_upload_successful_header, result['header'])
        self.assertEqual(strings.cancel_upload_successful_message, result['message'])
