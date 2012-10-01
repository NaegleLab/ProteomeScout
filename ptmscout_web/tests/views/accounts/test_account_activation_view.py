from mock import patch
from ptmscout.config import strings
from ptmscout.views.accounts.account_activation_view import \
    user_account_activation
from pyramid import testing
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockUser
import ptmscout.database.user as dbuser
import unittest

class UserAccountActivationTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
     
    def test_user_account_activation_should_fail_if_fields_not_specified(self):
        request = DummyRequest()
        
        request.GET['username'] = ""
        request.GET['token'] = ""
        
        result = user_account_activation(request)
        
        self.assertEqual(strings.account_activation_failed_header, result['header'])
        self.assertEqual(strings.account_activation_page_title, result['pageTitle'])
        self.assertEqual(request.application_url + "/register", result['redirect'])
        self.assertEqual(strings.account_activation_failed_message % (request.application_url + "/register"), result['message'])
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_account_activation_should_fail_if_username_does_not_exist(self, patch_getUser):
        patch_getUser.side_effect = dbuser.NoSuchUser()
        request = DummyRequest()
        
        request.GET['username'] = "username"
        request.GET['token'] = "token"
        
        result = user_account_activation(request)
        
        self.assertEqual(strings.account_activation_failed_header, result['header'])
        self.assertEqual(strings.account_activation_page_title, result['pageTitle'])
        self.assertEqual(request.application_url + "/register", result['redirect'])
        self.assertEqual(strings.account_activation_failed_message % (request.application_url + "/register"), result['message'])
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_account_activation_should_fail_if_user_token_incorrect(self, patch_getUser):
        patch_getUser.return_value = createMockUser("username", "email", "password", 0)
        request = DummyRequest()
        
        request.GET['username'] = "username"
        request.GET['token'] = "token"
        
        result = user_account_activation(request)
        
        self.assertEqual(strings.account_activation_failed_header, result['header'])
        self.assertEqual(strings.account_activation_page_title, result['pageTitle'])
        self.assertEqual(request.application_url + "/register", result['redirect'])
        self.assertEqual(strings.account_activation_failed_message % (request.application_url + "/register"), result['message'])
    
    @patch('ptmscout.database.user.User')
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_account_activation_should_succeed_if_all_parameters_correct(self, patch_getUser, patch_user):
        ptm_user = createMockUser("username", "email", "password", 0)
        patch_getUser.return_value = ptm_user
        request = DummyRequest()
        
        request.GET['username'] = "username"
        request.GET['token'] = ptm_user.activation_token
        
        result = user_account_activation(request)
        
        self.assertTrue(ptm_user.setActive.called, "User was not set as active")
        self.assertTrue(ptm_user.saveUser.called, "User was not saved")
        self.assertEqual(strings.account_activation_success_header, result['header'])
        self.assertEqual(strings.account_activation_page_title, result['pageTitle'])
        self.assertEqual(request.application_url + "/login", result['redirect'])
        self.assertEqual(strings.account_activation_success_message % (request.application_url+"/login"), result['message'])
