from pyramid import testing
from pyramid.httpexceptions import HTTPFound
from pyramid.testing import DummyRequest
import unittest
import urllib
from mock import patch, Mock
import ptmscout.database.user as dbuser
from ptmscout.user_management import user_login, user_logout, user_login_success,\
    user_registration_view, user_registration_success, user_account_activation,\
    forgot_password, process_forgot_password, manage_account, change_password,\
    change_password_success
import ptmscout.utils.crypto as crypto

class UserManagementTests(unittest.TestCase):
    TEST_USER_ID = 17
    
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
    def test_manage_account_should_display_account_info(self):
        request = DummyRequest()
        request.GET['reason'] = "reason"
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user        
        
        info = manage_account(request)
        self.assertEqual("Account Management", info['pageTitle'])
        self.assertEqual("username", info['username'])
        self.assertEqual("A User", info['fullname'])
        self.assertEqual("email", info['email'])
        self.assertEqual("institution", info['institution'])
        self.assertEqual("reason", info['reason'])
    
    @patch('ptmscout.utils.transactions.commit')
    def test_change_password_should_notify_on_success_and_update_db(self, patch_commit):
        request = DummyRequest()
        request.POST['old_pass'] = "password"
        request.POST['new_pass1'] = "newpassword"
        request.POST['new_pass2'] = "newpassword"
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        
        try:
            change_password(request)
            
            _, new_salted_password = crypto.saltedPassword("newpassword", ptm_user.salt)
            self.assertEqual(new_salted_password, ptm_user.salted_password)
            
            self.assertTrue(ptm_user.saveUser.called)
            self.assertTrue(patch_commit.called)
        except HTTPFound, f:
            self.assertEqual("http://example.com/change_password_success", f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
        
    def test_change_password_success_should_display_notification(self):
        request = DummyRequest()
        info = change_password_success(request)
        
        self.assertEqual("Password successfully changed.", info['message'])
        self.assertEqual("Success", info['header'])
        self.assertEqual("Change Password", info['pageTitle'])
        
    def test_change_password_should_fail_if_passwords_not_supplied(self):
        request = DummyRequest()
        request.POST['old_pass'] = ""
        request.POST['new_pass1'] = ""
        request.POST['new_pass2'] = ""
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        
        try:
            change_password(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/account?"+urllib.urlencode({'reason':"Form fields cannot be empty"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
    
    def test_change_password_should_fail_if_old_password_not_correct(self):
        request = DummyRequest()
        request.POST['old_pass'] = "wrongpass"
        request.POST['new_pass1'] = "newpassword"
        request.POST['new_pass2'] = "newpassword"
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        
        try:
            change_password(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/account?"+urllib.urlencode({'reason':"Supplied password was incorrect"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised") 
        
    def test_change_password_should_fail_if_passwords_do_not_match(self):
        request = DummyRequest()
        request.POST['old_pass'] = "password"
        request.POST['new_pass1'] = "newpass1"
        request.POST['new_pass2'] = "newpass2"
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        
        try:
            change_password(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/account?"+urllib.urlencode({'reason':"New password confirmation did not match"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")

    def test_forgot_password_should_display_forgot_password_form(self):
        request = DummyRequest()
        request.GET['email'] = "someemail@institute.edu"
        request.GET['reason'] = "reason"
        
        info = forgot_password(request)
        
        self.assertEqual("Forgotten Password Retrieval", info['pageTitle'])
        self.assertEqual("someemail@institute.edu", info['email'])
        self.assertEqual("reason", info['reason'])

    @patch('ptmscout.utils.transactions.commit')    
    @patch('ptmscout.database.user.getUserByEmail')
    @patch('ptmscout.utils.crypto.randomString')
    @patch('ptmscout.utils.mail.send_automail_message')
    def test_process_forgot_password_should_display_success_reset_password_and_send_email_when_user_exists(self, patch_sendMail, patch_crypto, patch_getUser, patch_commit):
        user_email = "a_valid_email@institute.edu"
        ptm_user = createUserForTest("username", user_email, "password", 1)
        patch_getUser.return_value = ptm_user
        new_password = "as24jdf945"
        patch_crypto.return_value = new_password
        
        request = DummyRequest()
        request.POST['email'] = user_email
        
        info = process_forgot_password(request)
        self.assertEqual("Forgotten Password Retrieval", info['pageTitle'])
        self.assertEqual("Password Reset Success", info['header'])
        self.assertEqual("Your username and a temporary password have been sent to your e-mail address", info['message'])
        
        self.assertTrue(patch_sendMail.called)
        
        login_url = request.application_url + "/login"
        account_url = request.application_url + "/account"
        
        password_reset_message = """A User,
        
        Your password in PTMScout has been reset, your new login credentials are:
        Username: %s
        Password: %s
        
        Please visit <a href="%s">PTMScout</a> to login.
        After logging in, your can change your password <a href="%s">here</a>.
        
        -PTMScout Administrator
        """ % ("username", new_password, login_url, account_url)
        
        patch_sendMail.assert_called_with(request, [user_email], "PTMScout password reset", password_reset_message)
        self.assertTrue(ptm_user.saveUser.called)
        self.assertTrue(patch_commit.called)
        
    @patch('ptmscout.database.user.getUserByEmail')
    def test_process_forgot_password_should_display_error_when_email_not_provided(self, patch_getUser):
        patch_getUser.side_effect = dbuser.NoSuchUser(username="a_username")
        request = DummyRequest()
        request.POST['email'] = ""
        
        try:
            process_forgot_password(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/forgot_password?"+urllib.urlencode({'email':"",'reason':"Form fields cannot be empty"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
            
    @patch('ptmscout.database.user.getUserByEmail')
    def test_process_forgot_password_should_display_error_when_no_such_user(self, patch_getUser):
        patch_getUser.side_effect = dbuser.NoSuchUser(username="a_username")
        request = DummyRequest()
        request.POST['email'] = "username@institute.edu"
        
        try:
            process_forgot_password(request)
        except HTTPFound, f:
            self.assertEqual("http://example.com/forgot_password?"+urllib.urlencode({'email':"username@institute.edu",'reason':"E-mail address does not match any user record"}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
                
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
    
    @patch('ptmscout.utils.mail.send_automail_message')
    @patch('ptmscout.utils.transactions.commit')
    @patch('ptmscout.database.user.User')
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_send_email_on_success_and_store_new_user(self, patch_getUser, patch_User, patch_commit, patch_automail):
        patch_getUser.side_effect = dbuser.NoSuchUser(username="a_username")
        
        request = DummyRequest()
        
        request.POST['username'] = "a_username"
        request.POST['pass1'] = "password1"
        request.POST['pass2'] = "password1"
        request.POST['email'] = "email@institute.edu"
        request.POST['name'] = "myname"
        request.POST['institution'] = "institute"
        
        result = user_registration_success(request)
        
        user_instance = patch_User.return_value
        self.assertTrue(patch_automail.called, "User was not notified of account activation by e-mail")
        self.assertTrue(user_instance.saveUser.called, "User was not saved to database")
        self.assertTrue(user_instance.createUser.called, "User object not initialized before DB commit")
        
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
        self.assertEqual("The specified account is not valid, please try <a href=\"http://example.com/register\">registering</a>", result['message'])
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_account_activation_should_fail_if_username_does_not_exist(self, patch_getUser):
        patch_getUser.side_effect = dbuser.NoSuchUser()
        request = DummyRequest()
        
        request.GET['username'] = "username"
        request.GET['token'] = "token"
        
        result = user_account_activation(request)
        
        self.assertEqual("Account Activation Failure", result['header'])
        self.assertEqual("Account Activation", result['pageTitle'])
        self.assertEqual("The specified account is not valid, please try <a href=\"http://example.com/register\">registering</a>", result['message'])
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_account_activation_should_fail_if_user_token_incorrect(self, patch_getUser):
        patch_getUser.return_value = createUserForTest("username", "email", "password", 0)
        request = DummyRequest()
        
        request.GET['username'] = "username"
        request.GET['token'] = "token"
        
        result = user_account_activation(request)
        
        self.assertEqual("Account Activation Failure", result['header'])
        self.assertEqual("Account Activation", result['pageTitle'])
        self.assertEqual("The specified account is not valid, please try <a href=\"http://example.com/register\">registering</a>", result['message'])
    
    @patch('ptmscout.utils.transactions.commit')
    @patch('ptmscout.database.user.User')
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_account_activation_should_succeed_if_all_parameters_correct(self, patch_getUser, patch_user, patch_commit):
        ptm_user = createUserForTest("username", "email", "password", 0)
        patch_getUser.return_value = ptm_user
        request = DummyRequest()
        
        request.GET['username'] = "username"
        request.GET['token'] = ptm_user.activation_token
        
        result = user_account_activation(request)
        
        self.assertTrue(ptm_user.setActive.called, "User was not set as active")
        self.assertTrue(ptm_user.saveUser.called, "User was not saved")
        self.assertTrue(patch_commit.called, "Database changed were not committed")
        self.assertEqual("Account Activation Succeeded", result['header'])
        self.assertEqual("Account Activation", result['pageTitle'])
        self.assertEqual("Your account is now active. Please <a href=\"http://example.com/login\">login</a>", result['message'])
    
    def test_salter(self):
        print crypto.saltedPassword("password", '41683edb7a')
        print crypto.saltedPassword("secret", '41683edb7a')
        
def createUserForTest(username, email, password, active):
    mock = Mock()
    mock.username = username
    mock.name = "A User"
    mock.email = email
    mock.institution = "institution"
    mock.salt, mock.salted_password = crypto.saltedPassword(password)  
    mock.activation_token = crypto.generateActivationToken()
    mock.id = UserManagementTests.TEST_USER_ID
    mock.active = active
    return mock
