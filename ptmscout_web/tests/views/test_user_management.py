from pyramid import testing
from pyramid.httpexceptions import HTTPFound
from pyramid.testing import DummyRequest
import unittest
import urllib
from mock import patch
import ptmscout.database.user as dbuser
from ptmscout.user_management import user_login, user_logout, user_login_success,\
    user_registration_view, user_registration_success, user_account_activation
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
        request.user = "something"
        value = user_logout(request)
        
        self.assertEqual(None, request.user)
        
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
        
        request.POST['username'] = "good_username"
        request.POST['password'] = "bad_password"

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
        
        request.POST['username'] = "good_username"
        request.POST['password'] = "good_password"

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
        patch_getUser.side_effect = dbuser.NoSuchUser(username="notauser")
        
        request = DummyRequest()
        
        request.POST['username'] = "notauser"
        request.POST['password'] = "password"
        
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
        
        request.POST['username'] = "good_username"
        request.POST['password'] = "good_password"

        value = user_login_success(request)
        
        self.assertEqual(patch_getUser.return_value, request.user)
        patch_security.assert_called_once_with(request, "good_username")
        self.assertEqual("Login", value['pageTitle'])
        self.assertEqual("Login Successful", value['header'])
        self.assertEqual("You have successfully logged in.", value['message'])
        
    def test_salter(self):
        print crypto.saltedPassword("password", '41683edb7a')
        
    def test_user_registration_view_should_display_correct_fields(self):
        request = DummyRequest()
        
        request.GET['username'] = "a_username"
        request.GET['email'] = "email@institute.edu"
        request.GET['name'] = "myname"
        request.GET['institution'] = "institute"
        request.GET['reason'] = "Passwords do not match"
        
        value = user_registration_view(request)
        
        self.assertEqual("a_username", value['username'])
        self.assertEqual("email@institute.edu", value['email'])
        self.assertEqual("myname", value['name'])
        self.assertEqual("institute", value['institution'])
        self.assertEqual("Passwords do not match", value['reason'])

    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_fail_and_redirect_when_empty_fields(self, patch_getUser):
        patch_getUser.side_effect = dbuser.NoSuchUser(username="a_username")
        request = DummyRequest()
        
        request.POST['username'] = "a_username"
        request.POST['pass1'] = "password1"
        request.POST['pass2'] = "password1"
        request.POST['email'] = "email@institute.edu"
        request.POST['name'] = "     "
        request.POST['institution'] = "institute"
        
        try:
            user_registration_success(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/register?"+urllib.urlencode({'username':"a_username", 'email': "email@institute.edu", 'name':"", 'institution':"institute",'reason':"Form fields cannot be empty"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")

    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_fail_and_redirect_when_email_is_not_valid(self, patch_getUser):
        patch_getUser.side_effect = dbuser.NoSuchUser(username="a_username")
        request = DummyRequest()
        
        request.POST['username'] = "a_username"
        request.POST['pass1'] = "password1"
        request.POST['pass2'] = "password1"
        request.POST['email'] = "invalidemail"
        request.POST['name'] = "myname"
        request.POST['institution'] = "institute"
        
        try:
            user_registration_success(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/register?"+urllib.urlencode({'username':"a_username", 'email': "invalidemail", 'name':"myname", 'institution':"institute",'reason':"Email address is invalid"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")

    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_fail_and_redirect_when_email_is_not_edu(self, patch_getUser):
        patch_getUser.side_effect = dbuser.NoSuchUser(username="a_username")
        request = DummyRequest()
        
        request.POST['username'] = "a_username"
        request.POST['pass1'] = "password1"
        request.POST['pass2'] = "password1"
        request.POST['email'] = "email@institute.com"
        request.POST['name'] = "myname"
        request.POST['institution'] = "institute"
        
        try:
            user_registration_success(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/register?"+urllib.urlencode({'username':"a_username", 'email': "email@institute.com", 'name':"myname", 'institution':"institute",'reason':"Email address must belong to .edu domain"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_fail_and_redirect_when_password_dont_match(self, patch_getUser):
        patch_getUser.side_effect = dbuser.NoSuchUser(username="a_username")
        request = DummyRequest()
        
        request.POST['username'] = "a_username"
        request.POST['pass1'] = "password1"
        request.POST['pass2'] = "password2"
        request.POST['email'] = "email@institute.edu"
        request.POST['name'] = "myname"
        request.POST['institution'] = "institute"
        
        try:
            user_registration_success(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/register?"+urllib.urlencode({'username':"a_username", 'email': "email@institute.edu", 'name':"myname", 'institution':"institute",'reason':"Password confirmation does not match"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
            
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_fail_and_redirect_when_password_too_short(self, patch_getUser):
        patch_getUser.side_effect = dbuser.NoSuchUser(username="a_username")            
        request = DummyRequest()
        
        request.POST['username'] = "a_username"
        request.POST['pass1'] = "P@$5"
        request.POST['pass2'] = "P@$5"
        request.POST['email'] = "email@institute.edu"
        request.POST['name'] = "myname"
        request.POST['institution'] = "institute"
        
        try:
            user_registration_success(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/register?"+urllib.urlencode({'username':"a_username", 'email': "email@institute.edu", 'name':"myname", 'institution':"institute",'reason':"Password must be at least 7 characters in length"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
        
    request = DummyRequest()
    @patch('ptmscout.database.user.getUserByUsername')    
    def test_user_registration_success_should_fail_when_user_already_exists(self, patch_getUser):
        patch_getUser.return_Value = createUserForTest("a_username", "email@institute.edu", "password", 1)
        
        request = DummyRequest()
        request.POST['username'] = "a_username"
        request.POST['pass1'] = "P@$5W0rd"
        request.POST['pass2'] = "P@$5W0rd"
        request.POST['email'] = "email@institute.edu"
        request.POST['name'] = "myname"
        request.POST['institution'] = "institute"
        
        try:
            user_registration_success(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/register?"+urllib.urlencode({'username':"a_username", 'email': "email@institute.edu", 'name':"myname", 'institution':"institute",'reason':"Username is already in use"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
        
    @patch('ptmscout.database.commit')
    @patch.object(dbuser.PTMUser, 'createUser')
    @patch.object(dbuser.PTMUser, 'saveUser')
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_send_email_on_success_and_store_new_user(self, patch_getUser, patch_saveUser, patch_createUser, patch_commit):
        patch_getUser.side_effect = dbuser.NoSuchUser(username="a_username")
        
        request = DummyRequest()
        
        request.POST['username'] = "a_username"
        request.POST['pass1'] = "password1"
        request.POST['pass2'] = "password1"
        request.POST['email'] = "email@institute.edu"
        request.POST['name'] = "myname"
        request.POST['institution'] = "institute"
        
        result = user_registration_success(request)

        self.assertTrue(patch_saveUser.called, "User was not saved to database")
        self.assertTrue(patch_createUser.called, "User object not initialized before DB commit")
        
        self.assertTrue(patch_commit.called, "Database changed were not committed")
        
        self.assertEqual("User Registration", result['pageTitle'])
        self.assertEqual("A confirmation e-mail has been sent to the specified e-mail address. Please check your e-mail to complete your registration.", result['message'])
        
        
    def test_user_account_activation_should_fail_if_fields_not_specified(self):
        request = DummyRequest()
        
        request.GET['username'] = ""
        request.GET['token'] = ""
        
        result = user_account_activation(request)
        
        self.assertEqual("Account Activation Failure", result['header'])
        self.assertEqual("Account Activation", result['pageTitle'])
        self.assertEqual("The specified account is not valid, please try <a href=\"http://www.example.com/register\">registering</a>", result['message'])
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_account_activation_should_fail_if_username_does_not_exist(self, patch_getUser):
        patch_getUser.side_effect = dbuser.NoSuchUser()
        request = DummyRequest()
        
        request.GET['username'] = "username"
        request.GET['token'] = "token"
        
        result = user_account_activation(request)
        
        self.assertEqual("Account Activation Failure", result['header'])
        self.assertEqual("Account Activation", result['pageTitle'])
        self.assertEqual("The specified account is not valid, please try <a href=\"http://www.example.com/register\">registering</a>", result['message'])
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_account_activation_should_fail_if_user_token_incorrect(self, patch_getUser):
        patch_getUser.return_value = createUserForTest("username", "email", "password", 0)
        request = DummyRequest()
        
        request.GET['username'] = "username"
        request.GET['token'] = "token"
        
        result = user_account_activation(request)
        
        self.assertEqual("Account Activation Failure", result['header'])
        self.assertEqual("Account Activation", result['pageTitle'])
        self.assertEqual("The specified account is not valid, please try <a href=\"http://www.example.com/register\">registering</a>", result['message'])
    
    @patch('ptmscout.database.commit')
    @patch.object(dbuser.PTMUser, 'saveUser')
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_account_activation_should_succeed_if_all_parameters_correct(self, patch_getUser, patch_saveUser, patch_commit):
        ptm_user = createUserForTest("username", "email", "password", 0)
        patch_getUser.return_value = ptm_user
        request = DummyRequest()
        
        request.GET['username'] = "username"
        request.GET['token'] = ptm_user.activation_token
        
        result = user_account_activation(request)
        
        self.assertTrue(patch_saveUser.called, "User was not saved")
        self.assertTrue(patch_commit.called, "Database changed were not committed")
        self.assertEqual(1, ptm_user.active)
        self.assertEqual("Account Activation Succeeded", result['header'])
        self.assertEqual("Account Activation", result['pageTitle'])
        self.assertEqual("Your account is now active. Please <a href=\"http://www.example.com/login\">login</a>", result['message'])
    
    
        
def createUserForTest(username, email, password, active):
    user = dbuser.PTMUser(username, "A User", email, "institution")
    user.createUser(password)
    user.id = UserManagementTests.TEST_USER_ID
    user.active = active
    return user
