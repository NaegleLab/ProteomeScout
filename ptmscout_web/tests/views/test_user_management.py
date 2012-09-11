from pyramid import testing
from pyramid.httpexceptions import HTTPFound
from pyramid.testing import DummyRequest
import unittest
import urllib
from mock import patch, Mock
from ptmscout.user_management import manage_account, change_password, change_password_success,\
    manage_experiments
import ptmscout.utils.crypto as crypto
from ptmscout import strings
from tests.views.mocking import createUserForTest, createMockExperiment,\
    createMockPermission

class UserManagementTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
    def test_my_experiments_should_show_experiments(self):
        request = DummyRequest()
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user

        exp1 = createMockExperiment(1, 0)
        exp2 = createMockExperiment(2, 0)
        exp3 = createMockExperiment(3, 0)
        
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        ptm_user.permissions.append(createMockPermission(ptm_user, exp2, 'owner'))
        ptm_user.permissions.append(createMockPermission(ptm_user, exp3, 'view'))
        
        info = manage_experiments(request)
        
        self.assertEqual([exp1, exp2], info['experiments'])
        self.assertEqual(strings.my_experiments_page_title, info['pageTitle'])
        
    def test_manage_account_should_display_account_info(self):
        request = DummyRequest()
        request.GET['reason'] = "reason"
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        
        info = manage_account(request)
        self.assertEqual(strings.account_management_page_title, info['pageTitle'])
        self.assertEqual("username", info['username'])
        self.assertEqual("A User", info['fullname'])
        self.assertEqual("email", info['email'])
        self.assertEqual("institution", info['institution'])
        self.assertEqual("reason", info['reason'])
    
    def test_change_password_should_notify_on_success_and_update_db(self):
        request = DummyRequest()
        request.POST['old_pass'] = "password"
        request.POST['new_pass1'] = "newpassword"
        request.POST['new_pass2'] = "newpassword"
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        
        reroute = change_password(request)
        
        _, new_salted_password = crypto.saltedPassword("newpassword", ptm_user.salt)
        self.assertEqual(new_salted_password, ptm_user.salted_password)
        
        self.assertTrue(ptm_user.saveUser.called)
        self.assertEqual(request.application_url + "/change_password_success", reroute.location)
        
    def test_change_password_success_should_display_notification(self):
        request = DummyRequest()
        info = change_password_success(request)
        
        self.assertEqual(strings.change_password_success_message, info['message'])
        self.assertEqual(strings.success_header, info['header'])
        self.assertEqual(strings.change_password_page_title, info['pageTitle'])
        
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
            self.assertEqual(request.application_url + "/account?"+urllib.urlencode({'reason':strings.failure_reason_form_fields_cannot_be_empty}), f.location)
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
            self.assertEqual(request.application_url + "/account?"+urllib.urlencode({'reason':strings.failure_reason_incorrect_password}), f.location)
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
            self.assertEqual(request.application_url + "/account?"+urllib.urlencode({'reason':strings.failure_reason_new_passwords_not_matching}), f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception HTTPFound, no exception raised")
