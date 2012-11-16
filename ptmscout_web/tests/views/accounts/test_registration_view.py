from mock import patch, Mock
from ptmscout.config import strings, settings
from ptmscout.database import user
from ptmscout.views.accounts.registration_view import user_registration_view, \
    validate_username, validate_password_length, validate_password_match, validate_email
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockUser
from ptmscout.utils import forms
from ptmscout.database.user import NoSuchUser
from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase

class UserRegistrationViewIntegrationTests(IntegrationTestCase):
    def test_view_should_allow_visitation(self):
        self.ptmscoutapp.get("/register", status=200)

class UserRegistrationViewTests(UnitTestCase):
    
    @patch('ptmscout.database.user.getUserByUsername')
    def test_validate_username_should_check_for_existing_user(self, patch_getUser):
        someschema = Mock(spec=forms.FormSchema)
        patch_getUser.return_value = createMockUser(username='someusername')
        
        result = validate_username('Username', 'someusername', someschema)
        
        self.assertEqual(strings.failure_reason_username_inuse, result)
        
    @patch('ptmscout.database.user.getUserByUsername')
    def test_validate_username_should_pass(self, patch_getUser):
        someschema = Mock(spec=forms.FormSchema)
        patch_getUser.side_effect = NoSuchUser()
        
        result = validate_username('Username', 'someusername', someschema)
        
        self.assertEqual(None, result)
    
    def test_validate_password_should_fail_when_too_short(self):
        result = validate_password_length('pass','2short', "someschema")
        self.assertEqual(strings.failure_reason_password_too_short % settings.MINIMUM_PASSWORD_LENGTH, result)

    def test_validate_password_match_should_fail_when_no_match(self):
        someschema = Mock(spec=forms.FormSchema)
        someschema.get_form_value.return_value = 'otherpass'
        
        result = validate_password_match('pass', 'password', someschema)
        
        someschema.get_form_value.assert_called_once_with('pass1')
        self.assertEqual(strings.failure_reason_new_passwords_not_matching, result)
    
    @patch('ptmscout.utils.mail.email_is_valid')
    def test_validate_email_should_fail_if_email_not_valid(self, patch_verify):
        patch_verify.return_value = False, False
        
        result = validate_email('email', 'notanemail', "someschema")
        
        patch_verify.assert_called_once_with('notanemail')
        self.assertEqual(strings.failure_reason_email_not_valid, result)
    
    @patch('ptmscout.utils.mail.email_is_valid')
    def test_validate_email_should_fail_if_email_not_correct_domain(self, patch_verify):
        patch_verify.return_value = True, False
        
        result = validate_email('email', 'notanemail@example.com', "someschema")
        
        patch_verify.assert_called_once_with('notanemail@example.com')
        self.assertEqual(strings.failure_reason_email_not_allowed, result)
    
    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.utils.forms.FormSchema')
    @patch('ptmscout.database.user.getUserByUsername')    
    def test_user_registration_success_should_fail_when_form_doesnt_validate(self, patch_getUser, patch_schema, patch_validator):
        patch_getUser.return_Value = createMockUser("a_username", "email@institute.edu", "password", 1)
        
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.POST['username'] = "a_username"
        request.POST['pass1'] = "P@$5W0rd"
        request.POST['pass2'] = "P@$5W0rd"
        request.POST['email'] = "email@institute.edu"
        request.POST['name'] = "myname"
        request.POST['institution'] = "institute"

        mock_validator = patch_validator.return_value
        mock_validator.validate.return_value = ["An error"]
        
        result = user_registration_view(request)
        
        self.assertEqual(["An error"], result['errors'])
        self.assertEqual(patch_schema.return_value, result['formrenderer'].schema)
        self.assertEqual(strings.user_registration_page_title, result['pageTitle'])
        
    @patch('ptmscout.utils.forms.FormValidator')
    @patch('ptmscout.utils.forms.FormSchema')
    @patch('ptmscout.utils.mail.send_automail_message')
    @patch('ptmscout.database.user.User')
    @patch('ptmscout.database.user.getUserByUsername')
    def test_user_registration_success_should_send_email_on_success_store_new_user_and_process_invitations(self, patch_getUser, patch_User, patch_automail, patch_schema, patch_validator):
        patch_getUser.side_effect = user.NoSuchUser(username="a_username")
        
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.POST['username'] = "a_username"
        request.POST['pass1'] = "password1"
        request.POST['pass2'] = "password1"
        request.POST['email'] = "email@institute.edu"
        request.POST['name'] = "myname"
        request.POST['institution'] = "institute"
        
        mock_validator = patch_validator.return_value
        mock_validator.validate.return_value = []
        
        f = user_registration_view(request)
        self.assertEqual(request.application_url + "/registration_success", f.location)
        
        user_instance = patch_User.return_value
        user_instance.processInvitations.assert_called_with()
        self.assertTrue(patch_automail.called, "User was not notified of account activation by e-mail")
        self.assertTrue(user_instance.saveUser.called, "User was not saved to database")
        self.assertTrue(user_instance.createUser.called, "User object not initialized before DB commit")
        
