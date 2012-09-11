from pyramid import testing
from pyramid.httpexceptions import HTTPFound, HTTPForbidden
from pyramid.testing import DummyRequest
import unittest
import urllib
from mock import patch, Mock
from ptmscout.user_management import manage_account, change_password, change_password_success,\
    manage_experiments, manage_experiment_permissions, publish_experiment
import ptmscout.utils.crypto as crypto
from ptmscout import strings
from tests.views.mocking import createUserForTest, createMockExperiment,\
    createMockPermission
from ptmscout.database.user import NoSuchUser

class UserManagementTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()

            
    def test_publish_experiment_should_throw_forbidden(self):
        request = DummyRequest()
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        request.matchdict['id'] = "100"
        
        try:
            publish_experiment(request)
        except HTTPForbidden:
            pass
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception HTTPForbidden was not thrown")
            
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_publish_experiment_should_show_success_already_public(self, patch_getExperiment):
        request = DummyRequest()
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        exp1.public = 1
        
        patch_getExperiment.return_value = exp1
        
        info = publish_experiment(request)
        
        self.assertEqual("true", info['confirm'])
        self.assertEqual(strings.publish_experiment_page_title, info['pageTitle'])
        self.assertEqual(strings.publish_experiment_already_message, info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(request.application_url + "/account/experiments", info['redirect'])
        
        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_publish_experiment_should_check_confirmation(self, patch_getExperiment):
        request = DummyRequest()
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        
        patch_getExperiment.return_value = exp1
        
        info = publish_experiment(request)
        
        self.assertEqual("false", info['confirm'])
        self.assertEqual(strings.publish_experiment_page_title, info['pageTitle'])
        self.assertEqual(strings.publish_experiment_confirm_message, info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(None, info['redirect'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_publish_experiment_should_change_public_flag(self, patch_getExperiment):    
        request = DummyRequest()
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        request.POST['confirm'] = "true"
        
        patch_getExperiment.return_value = exp1
        
        info = publish_experiment(request)
        
        exp1.makePublic.assert_called_once()
        exp1.saveExperiment.assert_called_once()
        
        self.assertEqual("true", info['confirm'])
        self.assertEqual(strings.publish_experiment_page_title, info['pageTitle'])
        self.assertEqual(strings.publish_experiment_success_message, info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(request.application_url + "/account/experiments", info['redirect'])
    
    @patch('ptmscout.database.user.getUserByEmail')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_manage_experiment_fails_to_add_user_if_no_such_user(self, patch_getExperiment, patch_getUser):
        request = DummyRequest()
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user

        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        
        patch_getExperiment.return_value = exp1
        
        user1 = createUserForTest("other", "email", "password", 1)
        
        exp1.permissions = [createMockPermission(ptm_user, exp1, 'owner'),
                            createMockPermission(user1, exp1, 'view')]
        
        
        request.POST['submitted'] = "1"
        request.POST['email'] = "emailaddr"
        
        patch_getUser.side_effect = NoSuchUser()
        
        info = manage_experiment_permissions(request)
        
        self.assertEqual([user1], info['users'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(strings.share_experiment_page_title, info['pageTitle'])
        self.assertEqual(strings.failure_reason_email_address_not_on_record, info['reason'])
    
    @patch('ptmscout.database.user.getUserByEmail')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_manage_experiment_adds_user_on_form_submission(self, patch_getExperiment, patch_getUser):
        request = DummyRequest()
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user

        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        
        patch_getExperiment.return_value = exp1
        
        user1 = createUserForTest("other", "email", "password", 1)
        
        exp1.permissions = [createMockPermission(ptm_user, exp1, 'owner'),
                            createMockPermission(user1, exp1, 'view')]
        
        
        request.POST['submitted'] = "1"
        request.POST['email'] = "emailaddr"
        new_user = createUserForTest("other2", "emailaddr", "password", 1)
        patch_getUser.return_value = new_user

        info = manage_experiment_permissions(request)
        
        self.assertEqual([user1, new_user], info['users'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(strings.share_experiment_page_title, info['pageTitle'])
        self.assertEqual(None, info['reason'])
        
        exp1.grantPermission.assert_called_once_with(new_user, 'view')
        exp1.saveExperiment.assert_called_once()
    
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_manage_experiment_shows_users(self, patch_getExperiment):
        request = DummyRequest()
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user

        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        
        patch_getExperiment.return_value = exp1
        
        user1 = createUserForTest("other", "email", "password", 1)
        
        exp1.permissions = [createMockPermission(ptm_user, exp1, 'owner'),
                            createMockPermission(user1, exp1, 'view')]
        
        info = manage_experiment_permissions(request)
        
        self.assertEqual([user1], info['users'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(strings.share_experiment_page_title, info['pageTitle'])
        self.assertEqual(None, info['reason'])

    
    def test_manage_experiment_permissions_throws_forbidden(self):
        request = DummyRequest()
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user

        exp1 = createMockExperiment(1, 0)
        exp2 = createMockExperiment(2, 0)
        exp3 = createMockExperiment(3, 0)
        
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        ptm_user.permissions.append(createMockPermission(ptm_user, exp2, 'owner'))
        ptm_user.permissions.append(createMockPermission(ptm_user, exp3, 'view'))
        
        request.matchdict['id'] = "%d" % (exp3.id)

        try:
            manage_experiment_permissions(request)
        except HTTPForbidden:
            pass
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception HTTPForbidden to be thrown")
        
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
        self.assertEqual(request.application_url + "/account", info['redirect'])
        
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
