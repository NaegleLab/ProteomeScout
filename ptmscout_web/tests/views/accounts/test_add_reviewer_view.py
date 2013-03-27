from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
from ptmscout.views.accounts import add_reviewer_view
from ptmscout.config import strings
from tests.views.mocking import createMockExperiment, createMockUser
from mock import patch, Mock

class TestAddReviewerView(UnitTestCase):
    
    def test_get_login_link(self):
        exp = createMockExperiment()
        request = DummyRequest()
        
        request.route_url = Mock()
        request.route_path = Mock()
        
        request.route_url.return_value = 'http://example.com/login'
        request.route_path.return_value = '/experiments/26'
        
        result = add_reviewer_view.get_login_link(exp, request)
        
        self.assertEqual('http://example.com/login?redirect=%2Fexperiments%2F26', result)
        
        request.route_url.assert_called_once_with('login')
        request.route_path.assert_called_once_with('experiment', id=exp.id)
    
    @patch('ptmscout.views.accounts.add_reviewer_view.get_login_link')
    @patch('ptmscout.views.accounts.add_reviewer_view.create_review_user')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_add_reviewer_view_should_display_auth_token_link_after_confirm(self, patch_getExp, patch_createUser, patch_getLink):
        request = DummyRequest()
        request.matchdict['id'] = '34'
        request.POST['confirm'] = "true"
        request.user = createMockUser()
        
        exp = createMockExperiment()
        patch_getExp.return_value = exp
        patch_createUser.return_value = 'username', 'password'
        patch_getLink.return_value = 'login_link'
        
        result = add_reviewer_view.create_reviewer_account(request)
        
        patch_createUser.assert_called_once_with(exp)
        patch_getExp.assert_called_once_with(34, request.user)
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.create_reviewer_account_created_header, result['header'])
        self.assertEqual(strings.create_reviewer_account_created_message % ('username', 'password', 'login_link'), result['message'])
        self.assertEqual(strings.create_reviewer_account_page_title, result['pageTitle'])
        self.assertEqual(None, result['redirect'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_add_reviewer_view_should_display_confirmation_message(self, patch_getExp):
        request = DummyRequest()
        request.matchdict['id'] = '34'
        request.user = createMockUser()
        
        exp = createMockExperiment()
        patch_getExp.return_value = exp

        result = add_reviewer_view.create_reviewer_account(request)
        
        patch_getExp.assert_called_once_with(34, request.user)
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(strings.create_reviewer_account_confirm_header, result['header'])
        self.assertEqual(strings.create_reviewer_account_confirm_message, result['message'])
        self.assertEqual(strings.create_reviewer_account_page_title, result['pageTitle'])
        self.assertEqual(None, result['redirect'])

        
class IntegrationTestAddReviewerView(IntegrationTestCase):
    def test_view_integration(self):
        self.bot.acquire_experiments([26])
        result = self.ptmscoutapp.get("/account/experiments/26/add_reviewer", status=200)
        result = result.forms[0].submit()
        
        result.mustcontain('username: reviewer_')
        result.mustcontain('Anonymous reviewers may login to view your data by visiting the following address in their browser: http://localhost/login?redirect=%2Fexperiments%2F26')