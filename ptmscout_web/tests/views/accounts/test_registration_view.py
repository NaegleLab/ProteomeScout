from mock import patch
from ptmscout.config import strings
from ptmscout.database import user
from ptmscout.views.accounts.registration_view import user_registration_view, \
    user_registration_success
from pyramid import testing
from pyramid.httpexceptions import HTTPFound
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockUser
import unittest
import urllib

class UserRegistrationViewTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
    def test_user_registration_view_should_display_correct_fields(self):
        request = DummyRequest()
        
        request.GET['username'] = "a_username"
        request.GET['email'] = "email@institute.edu"
        request.GET['name'] = "myname"
        request.GET['institution'] = "institute"
        request.GET['reason'] = "reason"
        
        value = user_registration_view(request)
        
        self.assertEqual("a_username", value['username'])
        self.assertEqual("email@institute.edu", value['email'])
        self.assertEqual("myname", value['name'])
        self.assertEqual("institute", value['institution'])
        self.assertEqual("reason", value['reason'])

    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_fail_and_redirect_when_empty_fields(self, patch_getUser):
        patch_getUser.side_effect = user.NoSuchUser(username="a_username")
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
            self.assertEqual(request.application_url + "/register?"+urllib.urlencode({'username':"a_username", 'email': "email@institute.edu", 'name':"", 'institution':"institute",'reason':strings.failure_reason_form_fields_cannot_be_empty}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")

    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_fail_and_redirect_when_email_is_not_valid(self, patch_getUser):
        patch_getUser.side_effect = user.NoSuchUser(username="a_username")
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
            self.assertEqual(request.application_url + "/register?"+urllib.urlencode({'username':"a_username", 'email': "invalidemail", 'name':"myname", 'institution':"institute",'reason':strings.failure_reason_email_not_valid}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")

    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_fail_and_redirect_when_email_is_not_edu(self, patch_getUser):
        patch_getUser.side_effect = user.NoSuchUser(username="a_username")
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
            self.assertEqual(request.application_url + "/register?"+urllib.urlencode({'username':"a_username", 'email': "email@institute.com", 'name':"myname", 'institution':"institute",'reason':strings.failure_reason_email_not_academic}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_fail_and_redirect_when_password_dont_match(self, patch_getUser):
        patch_getUser.side_effect = user.NoSuchUser(username="a_username")
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
            self.assertEqual(request.application_url + "/register?"+urllib.urlencode({'username':"a_username", 'email': "email@institute.edu", 'name':"myname", 'institution':"institute",'reason':strings.failure_reason_new_passwords_not_matching}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
            
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_fail_and_redirect_when_password_too_short(self, patch_getUser):
        patch_getUser.side_effect = user.NoSuchUser(username="a_username")            
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
            self.assertEqual(request.application_url + "/register?"+urllib.urlencode({'username':"a_username", 'email': "email@institute.edu", 'name':"myname", 'institution':"institute",'reason':strings.failure_reason_password_too_short % (7)}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
        
    request = DummyRequest()
    @patch('ptmscout.database.user.getUserByUsername')    
    def test_user_registration_success_should_fail_when_user_already_exists(self, patch_getUser):
        patch_getUser.return_Value = createMockUser("a_username", "email@institute.edu", "password", 1)
        
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
            self.assertEqual(request.application_url + "/register?"+urllib.urlencode({'username':"a_username", 'email': "email@institute.edu", 'name':"myname", 'institution':"institute",'reason':strings.failure_reason_username_inuse}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
    
    @patch('ptmscout.utils.mail.send_automail_message')
    @patch('ptmscout.database.user.User')
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_send_email_on_success_store_new_user_and_process_invitations(self, patch_getUser, patch_User, patch_automail):
        patch_getUser.side_effect = user.NoSuchUser(username="a_username")
        
        request = DummyRequest()
        
        request.POST['username'] = "a_username"
        request.POST['pass1'] = "password1"
        request.POST['pass2'] = "password1"
        request.POST['email'] = "email@institute.edu"
        request.POST['name'] = "myname"
        request.POST['institution'] = "institute"
        
        result = user_registration_success(request)
        
        user_instance = patch_User.return_value
        user_instance.processInvitations.assert_called_with()
        self.assertTrue(patch_automail.called, "User was not notified of account activation by e-mail")
        self.assertTrue(user_instance.saveUser.called, "User was not saved to database")
        self.assertTrue(user_instance.createUser.called, "User object not initialized before DB commit")
        
        self.assertEqual(strings.user_registration_page_title, result['pageTitle'])
        self.assertEqual(strings.user_registration_success_header, result['header'])
        self.assertEqual(strings.user_registration_success_message, result['message'])
        self.assertEqual(request.application_url + "/login", result['redirect'])
        