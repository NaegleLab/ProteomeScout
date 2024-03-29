from mock import patch
from ptmscout.config import strings
from ptmscout.views.accounts.login_view import user_login, user_logout, \
    user_login_success
from pyramid.httpexceptions import HTTPFound
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockUser
import urllib
from ptmscout.database import user
from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase
import datetime

class UserLoginViewIntegration(IntegrationTestCase):
    def test_login_with_redirect(self):
        self.bot.logout()
        
        result = self.ptmscoutapp.get('/login?origin=%2Fexperiments%2F26')
        result.form.set( 'username', self.bot.username )
        result.form.set( 'password', 'nottherightpassword' )
        
        result = result.form.submit().follow()
        result.form.set( 'username', self.bot.username )
        result.form.set( 'password', self.bot.password )
        
        result = result.form.submit()
#        result.showbrowser()
        result.mustcontain('<meta http-equiv="refresh" content="3; url=http://localhost/experiments/26" />')

    def test_login_expired_account(self):
        self.bot.user.expiration = datetime.datetime.now()
        self.bot.user.saveUser()
        self.bot.logout()
        
        result = self.ptmscoutapp.get('/login?origin=%2Fexperiments%2F26')
        result.form.set( 'username', self.bot.username )
        result.form.set( 'password', self.bot.password )
        
        result = result.form.submit().follow()
        
        result.mustcontain(strings.failure_reason_account_expired)

class UserLoginViewTests(UnitTestCase):
    def test_login_should_display_login_page(self):
        
        request = DummyRequest()
    
        info = user_login(request)
        
        self.assertEqual(strings.login_page_title, info['pageTitle'])
        self.assertEqual(None, info['reason'])
        
    @patch('pyramid.security.forget')
    def test_user_logout_should_display_logout_message(self, patch_security):
        request = DummyRequest()
        request.user = "something"
        value = user_logout(request)
        
        self.assertEqual(None, request.user)
        
        patch_security.assert_called_once_with(request)
        self.assertEqual(strings.logout_page_title, value['pageTitle'])
        self.assertEqual(strings.logout_page_header, value['header'])
        self.assertEqual(strings.logout_page_message, value['message'])
        self.assertEqual(request.application_url, value['redirect'])
        
    def test_login_should_display_error_with_reason_and_populate_username(self):
        request = DummyRequest()
        
        request.GET['username'] = "a_username"
        request.GET['reason'] = "you didn't do it right"
        
        info = user_login(request)
        
        self.assertEqual(strings.login_page_title, info['pageTitle'])
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
            self.assertEqual(request.application_url + "/login?"+urllib.urlencode({'username':"a_username",'reason':strings.failure_reason_form_fields_cannot_be_empty}), f.location)
    

    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_process_login_should_fail_if_credentials_incorrect(self, patch_getUser):
        patch_getUser.return_value = createMockUser(username="good_username", email="user@institute.edu", password="good_password", active=1)
        
        request = DummyRequest()
        
        request.POST['username'] = "good_username"
        request.POST['password'] = "bad_password"

        try:
            user_login_success(request)
        except HTTPFound, f:
            self.assertEqual(request.application_url + "/login?"+urllib.urlencode({'username':"good_username",'reason':strings.failure_reason_incorrect_credentials}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_process_login_should_fail_if_account_inactive(self, patch_getUser):
        patch_getUser.return_value = createMockUser("good_username", "user@institute.edu", "good_password", 0)
        
        request = DummyRequest()
        
        request.POST['username'] = "good_username"
        request.POST['password'] = "good_password"

        try:
            user_login_success(request)
        except HTTPFound, f:
            self.assertEqual(request.application_url + "/login?"+urllib.urlencode({'username':"good_username",'reason':strings.failure_reason_inactive_account}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_process_login_should_fail_if_no_such_user(self, patch_getUser):
        patch_getUser.side_effect = user.NoSuchUser(username="notauser")
        
        request = DummyRequest()
        
        request.POST['username'] = "notauser"
        request.POST['password'] = "password"
        
        try:
            user_login_success(request)
        except HTTPFound, f:
            self.assertEqual(request.application_url + "/login?"+urllib.urlencode({'username':"notauser",'reason':strings.failure_reason_incorrect_credentials}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
    
    @patch('pyramid.security.remember')
    @patch('ptmscout.database.user.getUserByUsername')        
    def test_process_login_should_succeed_if_credentials_correct(self, patch_getUser, patch_security):
        patch_getUser.return_value = createMockUser("good_username", "user@institute.edu", "good_password", 1)
        
        request = DummyRequest()
        
        request.POST['username'] = "good_username"
        request.POST['password'] = "good_password"

        value = user_login_success(request)
        
        self.assertEqual(patch_getUser.return_value, request.user)
        patch_security.assert_called_once_with(request, "good_username")
        self.assertEqual(strings.login_page_title, value['pageTitle'])
        self.assertEqual(strings.login_page_success_header, value['header'])
        self.assertEqual(strings.login_page_success_message, value['message'])
        self.assertEqual(request.application_url, value['redirect'])
        
