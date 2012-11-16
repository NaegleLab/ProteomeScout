from mock import patch
from ptmscout.config import strings
from ptmscout.database import user
from ptmscout.views.accounts.forgot_password_view import forgot_password, \
    process_forgot_password
from pyramid.httpexceptions import HTTPFound
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockUser
import urllib
from tests.PTMScoutTestCase import UnitTestCase


class UserForgotPasswordViewTests(UnitTestCase):
   
    def test_forgot_password_should_display_forgot_password_form(self):
        request = DummyRequest()
        request.GET['email'] = "someemail@institute.edu"
        request.GET['reason'] = "reason"
        
        info = forgot_password(request)
        
        self.assertEqual(strings.forgotten_password_page_title, info['pageTitle'])
        self.assertEqual("someemail@institute.edu", info['email'])
        self.assertEqual("reason", info['reason'])

    @patch('ptmscout.database.user.getUserByEmail')
    @patch('ptmscout.utils.crypto.randomString')
    @patch('ptmscout.utils.mail.send_automail_message')
    def test_process_forgot_password_should_display_success_reset_password_and_send_email_when_user_exists(self, patch_sendMail, patch_crypto, patch_getUser):
        user_email = "a_valid_email@institute.edu"
        ptm_user = createMockUser(email=user_email, password="some_pass")
        patch_getUser.return_value = ptm_user
        new_password = "as24jdf945"
        patch_crypto.return_value = new_password
        
        request = DummyRequest()
        request.POST['email'] = user_email
        
        info = process_forgot_password(request)
        self.assertEqual(strings.forgotten_password_page_title, info['pageTitle'])
        self.assertEqual(strings.forgotten_password_success_header, info['header'])
        self.assertEqual(strings.forgotten_password_success_message, info['message'])
        self.assertEqual(request.application_url + "/login", info['redirect'])
        
        self.assertTrue(patch_sendMail.called)
        
        login_url = request.application_url + "/login"
        account_url = request.application_url + "/account"
        
        password_reset_message = strings.forgotten_password_email_message % (ptm_user.name, ptm_user.username, new_password, login_url, account_url)
        
        patch_sendMail.assert_called_with(request, [user_email], 
                                          strings.forgotten_password_email_subject, password_reset_message)
        self.assertTrue(ptm_user.saveUser.called)
        
    @patch('ptmscout.database.user.getUserByEmail')
    def test_process_forgot_password_should_display_error_when_email_not_provided(self, patch_getUser):
        patch_getUser.side_effect = user.NoSuchUser(username="a_username")
        request = DummyRequest()
        request.POST['email'] = ""
        
        try:
            process_forgot_password(request)
        except HTTPFound, f:
            self.assertEqual(request.application_url + "/forgot_password?"+urllib.urlencode({'email':"",'reason':strings.failure_reason_form_fields_cannot_be_empty}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
            
    @patch('ptmscout.database.user.getUserByEmail')
    def test_process_forgot_password_should_display_error_when_no_such_user(self, patch_getUser):
        patch_getUser.side_effect = user.NoSuchUser(username="a_username")
        request = DummyRequest()
        request.POST['email'] = "username@institute.edu"
        
        try:
            process_forgot_password(request)
        except HTTPFound, f:
            self.assertEqual(request.application_url + "/forgot_password?"+urllib.urlencode({'email':"username@institute.edu",'reason':strings.failure_reason_email_address_not_on_record}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
                