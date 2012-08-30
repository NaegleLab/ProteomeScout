from pyramid import testing
from pyramid.httpexceptions import HTTPFound
from pyramid.testing import DummyRequest
import unittest
import urllib
from mock import patch
from ptmscout.database.user import NoSuchUser, PTMUser
from ptmscout.user_management import user_login, user_logout, user_login_success
import ptmscout.utils.crypto as crypto

class UserManagementTests(unittest.TestCase):
    TEST_USER_ID = 17
    
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()

    def test_login_should_display_login_page(self):
        
        request = DummyRequest()
    
        info = user_login(request)
        
        self.assertEqual('Login', info['pageTitle'])
        self.assertEqual(None, info['reason'])
        
    @patch('pyramid.security.forget')
    def test_user_logout_should_display_logout_message(self, patch_security):
        request = DummyRequest()
        value = user_logout(request)
        
        patch_security.assert_called_once_with(request)
        self.assertEqual("Logout", value['pageTitle'])
        self.assertEqual("Logout Successful", value['header'])
        self.assertEqual("You have successfully logged out.", value['message'])
        
    def test_login_should_display_error_with_reason_and_populate_username(self):
        request = DummyRequest()
        
        request.GET['username'] = "a_username"
        request.GET['reason'] = "you didn't do it right"
        
        info = user_login(request)
        
        self.assertEqual('Login', info['pageTitle'])
        self.assertEqual('a_username', info['username'])
        self.assertEqual("you didn't do it right", info['reason'])
        
    def test_process_login_should_redirect_to_login_with_error_on_fields_missing(self):
        request = DummyRequest()
        
        request.POST['username'] = "a_username"
        request.POST['password'] = ""
        
        try:
            user_login_success(request)
            self.fail("Expected exception HTTPFound, no exception raised")
        except HTTPFound, f:
            self.assertEqual("http://example.com/login?"+urllib.urlencode({'username':"a_username",'reason':"All fields are required"}), f.location)
    

    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_process_login_should_fail_if_credentials_incorrect(self, patch_getUser):
        patch_getUser.return_value = createUserForTest("good_username", "user@institute.edu", "good_password", 1)
        
        request = DummyRequest()
        
        request.GET['username'] = "good_username"
        request.GET['password'] = "bad_password"

        try:
            user_login_success(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/login?"+urllib.urlencode({'username':"good_username",'reason':"Credentials incorrect"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_process_login_should_fail_if_account_inactive(self, patch_getUser):
        patch_getUser.return_value = createUserForTest("good_username", "user@institute.edu", "good_password", 0)
        
        request = DummyRequest()
        
        request.GET['username'] = "good_username"
        request.GET['password'] = "good_password"

        try:
            user_login_success(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/login?"+urllib.urlencode({'username':"good_username",'reason':"Account has not been activated"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_process_login_should_fail_if_no_such_user(self, patch_getUser):
        patch_getUser.side_effect = NoSuchUser(username="notauser")
        
        request = DummyRequest()
        
        request.GET['username'] = "notauser"
        request.GET['password'] = "password"
        
        try:
            user_login_success(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/login?"+urllib.urlencode({'username':"notauser",'reason':"Credentials incorrect"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
    
    @patch('pyramid.security.remember')
    @patch('ptmscout.database.user.getUserByUsername')        
    def test_process_login_should_succeed_if_credentials_correct(self, patch_getUser, patch_security):
        patch_getUser.return_value = createUserForTest("good_username", "user@institute.edu", "good_password", 1)
        
        request = DummyRequest()
        
        request.GET['username'] = "good_username"
        request.GET['password'] = "good_password"

        value = user_login_success(request)
        
        patch_security.assert_called_once_with(request, UserManagementTests.TEST_USER_ID)
        self.assertEqual("Login", value['pageTitle'])
        self.assertEqual("Login Successful", value['header'])
        self.assertEqual("You have successfully logged in.", value['message'])
        
    def test_salter(self):
        print crypto.saltedPassword("password", '41683edb7a')
        
def createUserForTest(username, email, password, active):
    user = PTMUser(username, "A User", email, "institution")
    user.createUser(password)
    user.id = UserManagementTests.TEST_USER_ID
    user.active = active
    return user
