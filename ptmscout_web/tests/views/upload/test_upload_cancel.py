from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.views.upload.upload_cancel import cancel_upload_view
from tests.views.mocking import createMockSession, createMockUser,\
    createMockExperiment
from pyramid.testing import DummyRequest
from ptmscout.config import strings
from mock import patch

class TestCancelView(UnitTestCase):


    @patch('ptmscout.database.upload.getExperimentById')
    @patch('ptmscout.database.upload.getSessionById')
    def test_view(self, patch_session, patch_experiment):
        user = createMockUser()
        exp = createMockExperiment()
        session = createMockSession(user,experiment_id=exp.id)
        
        patch_session.return_value = session  
        patch_experiment.return_value = exp
        
        request = DummyRequest()
        request.matchdict['id'] = '247'
        request.user = user
        
        result = cancel_upload_view(request)

        patch_session.assert_called_once_with(247, user)
        patch_experiment.assert_called_once_with(30)
        session.delete.assert_called_once_with()
        
        self.assertEqual(strings.cancel_upload_successful_page_title, result['pageTitle'])
        self.assertEqual(strings.cancel_upload_successful_header, result['header'])
        self.assertEqual(strings.cancel_upload_successful_message, result['message'])


    @patch('ptmscout.database.upload.getSessionById')
    def test_view_should_delete_session(self, patch_session):
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
