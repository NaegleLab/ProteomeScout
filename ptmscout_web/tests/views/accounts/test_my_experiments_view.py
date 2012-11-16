from mock import patch
from ptmscout.config import strings
from ptmscout.database.user import NoSuchUser
from ptmscout.views.accounts.my_experiments_view import confirm_invite_user, \
    privatize_experiment, publish_experiment, manage_experiment_permissions, \
    manage_experiments
from pyramid.httpexceptions import HTTPFound
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockUser, createMockExperiment, \
    createMockPermission
from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase

class MyExperimentsViewIntegrationTests(IntegrationTestCase):
    
    def test_integration(self):
        from ptmscout.database import experiment
        
        exp = experiment.getExperimentById(26, None, False)
        exp.status = 'loading'
        exp.grantPermission(self.bot.user, 'owner')
        exp.saveExperiment()
        
        self.bot.login()
        
        result = self.ptmscoutapp.get('/account/experiments')
        
        result.mustcontain(exp.name)
        result.mustcontain('Status')
        result.mustcontain('loading')
        

class MyExperimentsViewTests(UnitTestCase):
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_confirm_invite_should_fail_for_bad_domain(self, patch_getExperiment):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 1)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        invited_email = "inviteduser@institute.com"
        request.GET['email'] = invited_email
                
        patch_getExperiment.return_value = exp1
        
        info = confirm_invite_user(request)
        
        self.assertEqual(True, info['confirm'])
        self.assertEqual("Failure", info['header'])
        self.assertEqual(strings.user_invite_page_title, info['pageTitle'])
        self.assertEqual(strings.failure_reason_email_not_allowed, info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(request.application_url + "/account/experiments", info['redirect'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_confirm_invite_should_fail_for_bad_email(self, patch_getExperiment):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 1)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        invited_email = "inviteduser@institute"
        request.GET['email'] = invited_email
                
        patch_getExperiment.return_value = exp1
        
        info = confirm_invite_user(request)
        
        self.assertEqual(True, info['confirm'])
        self.assertEqual("Failure", info['header'])
        self.assertEqual(strings.user_invite_page_title, info['pageTitle'])
        self.assertEqual(strings.failure_reason_email_not_valid, info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(request.application_url + "/account/experiments", info['redirect'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_confirm_invite_should_check_confirmation(self, patch_getExperiment):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 1)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        invited_email = "inviteduser@institute.edu"
        request.GET['email'] = invited_email
                
        patch_getExperiment.return_value = exp1
        
        info = confirm_invite_user(request)
        
        self.assertEqual(False, info['confirm'])
        self.assertEqual("Confirm", info['header'])
        self.assertEqual(strings.user_invite_page_title, info['pageTitle'])
        self.assertEqual(strings.user_invite_confirm % (invited_email), info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(None, info['redirect'])

    @patch('ptmscout.database.permissions.Invitation')
    @patch('ptmscout.utils.mail.send_automail_message')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_confirm_invite_should_send_invitation(self, patch_getExperiment, patch_mail, patch_invite):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 1)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        request.POST['confirm'] = "true"
        
        invited_email = "inviteduser@institute.edu"
        request.GET['email'] = invited_email
                
        patch_getExperiment.return_value = exp1
        
        info = confirm_invite_user(request)
        
        inst = patch_invite.return_value
        inst.saveInvitation.assert_called_with()
        
        expected_email_body = strings.user_invite_email_message % (invited_email, ptm_user.name, exp1.name, request.application_url + "/register?email=" + invited_email)
        expected_email_body = expected_email_body.replace("\n", "\\n").replace("'","\\'")
        
        self.assertIn(strings.user_invite_email_subject % (ptm_user.name), str(patch_mail.call_args))
        self.assertIn(expected_email_body, str(patch_mail.call_args))
        self.assertEqual(True, info['confirm'])
        self.assertEqual("Success", info['header'])
        self.assertEqual(strings.user_invite_page_title, info['pageTitle'])
        self.assertEqual(strings.user_invited % (invited_email), info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(request.application_url + "/account/experiments", info['redirect'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_privatize_experiment_should_show_success_already_private(self, patch_getExperiment):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        exp1.public = 0
        
        patch_getExperiment.return_value = exp1
        
        info = privatize_experiment(request)
        
        self.assertEqual(True, info['confirm'])
        self.assertEqual("Success", info['header'])
        self.assertEqual(strings.privatize_experiment_page_title, info['pageTitle'])
        self.assertEqual(strings.privatize_experiment_already_message, info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(request.application_url + "/account/experiments", info['redirect'])
        
        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_privatize_experiment_should_check_confirmation(self, patch_getExperiment):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 1)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        
        patch_getExperiment.return_value = exp1
        
        info = privatize_experiment(request)
        
        self.assertEqual(False, info['confirm'])
        self.assertEqual("Confirm", info['header'])
        self.assertEqual(strings.privatize_experiment_page_title, info['pageTitle'])
        self.assertEqual(strings.privatize_experiment_confirm_message, info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(None, info['redirect'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_privatize_experiment_should_change_public_flag(self, patch_getExperiment):    
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 1)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        request.POST['confirm'] = "true"
        
        patch_getExperiment.return_value = exp1
        
        info = privatize_experiment(request)
        
        self.assert_(exp1.makePrivate.called)
        self.assert_(exp1.saveExperiment.called)
        
        self.assertEqual(True, info['confirm'])
        self.assertEqual("Success", info['header'])
        self.assertEqual(strings.privatize_experiment_page_title, info['pageTitle'])
        self.assertEqual(strings.privatize_experiment_success_message, info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(request.application_url + "/account/experiments", info['redirect'])
    

    @patch('ptmscout.database.experiment.getExperimentById')
    def test_publish_experiment_should_show_success_already_public(self, patch_getExperiment):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 1)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        
        patch_getExperiment.return_value = exp1
        
        info = publish_experiment(request)
        
        self.assertEqual(True, info['confirm'])
        self.assertEqual("Success", info['header'])
        self.assertEqual(strings.publish_experiment_page_title, info['pageTitle'])
        self.assertEqual(strings.publish_experiment_already_message, info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(request.application_url + "/account/experiments", info['redirect'])
        
        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_publish_experiment_should_check_confirmation(self, patch_getExperiment):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        
        patch_getExperiment.return_value = exp1
        
        info = publish_experiment(request)
        
        self.assertEqual(False, info['confirm'])
        self.assertEqual("Confirm", info['header'])
        self.assertEqual(strings.publish_experiment_page_title, info['pageTitle'])
        self.assertEqual(strings.publish_experiment_confirm_message, info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(None, info['redirect'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_publish_experiment_should_change_public_flag(self, patch_getExperiment):    
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        
        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        request.POST['confirm'] = "true"
        
        patch_getExperiment.return_value = exp1
        
        info = publish_experiment(request)
        
        exp1.makePublic.assert_called_once()
        exp1.saveExperiment.assert_called_once()
        
        self.assertEqual(True, info['confirm'])
        self.assertEqual("Success", info['header'])
        self.assertEqual(strings.publish_experiment_page_title, info['pageTitle'])
        self.assertEqual(strings.publish_experiment_success_message, info['message'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(request.application_url + "/account/experiments", info['redirect'])
    
    @patch('ptmscout.database.user.getUserByEmail')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_manage_experiment_redirects_to_invite_if_no_such_user(self, patch_getExperiment, patch_getUser):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user

        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        
        patch_getExperiment.return_value = exp1
        
        user1 = createMockUser("other", "email", "password", 1)
        
        exp1.permissions = [createMockPermission(ptm_user, exp1, 'owner'),
                            createMockPermission(user1, exp1, 'view')]
        
        request.POST['submitted'] = "1"
        request.POST['email'] = "emailaddr"
        
        patch_getUser.side_effect = NoSuchUser()

        try:
            manage_experiment_permissions(request)
        except HTTPFound, f:
            self.assertEqual(request.application_url + "/account/experiments/" + str(exp1.id) + "/invite?email=emailaddr", f.location)
        except Exception, e:
            self.fail("Unexpected exception thrown: " + str(e))
        else:
            self.fail("Expected exception: HTTPFound")
    
    
    @patch('ptmscout.utils.mail.send_automail_message')
    @patch('ptmscout.database.user.getUserByEmail')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_manage_experiment_adds_user_on_form_submission(self, patch_getExperiment, patch_getUser, patch_sendMail):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user

        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        
        patch_getExperiment.return_value = exp1
        
        user1 = createMockUser("other", "email", "password", 1)
        
        exp1.permissions = [createMockPermission(ptm_user, exp1, 'owner'),
                            createMockPermission(user1, exp1, 'view')]
        
        
        request.POST['submitted'] = "1"
        request.POST['email'] = "emailaddr"
        new_user = createMockUser("other2", "emailaddr", "password", 1)
        patch_getUser.return_value = new_user

        info = manage_experiment_permissions(request)
        
        patch_mail_result = str(patch_sendMail.call_args).replace("\\n", "\n").replace("\\'", "'")
        self.assertIn(strings.user_invite_email_subject % ptm_user.name, patch_mail_result)
        self.assertIn(strings.user_invite_email_message % (new_user.name, ptm_user.name, exp1.name, request.application_url + "/login"), patch_mail_result)
        self.assertTrue(patch_sendMail.called)
        
        self.assertEqual([user1, new_user], info['users'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(strings.share_experiment_page_title, info['pageTitle'])
        self.assertEqual(None, info['reason'])
        
        exp1.grantPermission.assert_called_once_with(new_user, 'view')
        exp1.saveExperiment.assert_called_once()
    
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_manage_experiment_shows_users(self, patch_getExperiment):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user

        exp1 = createMockExperiment(1, 0)
        ptm_user.permissions.append(createMockPermission(ptm_user, exp1, 'owner'))
        request.matchdict['id'] = "%d" % (exp1.id)
        
        patch_getExperiment.return_value = exp1
        
        user1 = createMockUser("other", "email", "password", 1)
        
        exp1.permissions = [createMockPermission(ptm_user, exp1, 'owner'),
                            createMockPermission(user1, exp1, 'view')]
        
        info = manage_experiment_permissions(request)
        
        self.assertEqual([user1], info['users'])
        self.assertEqual(exp1, info['experiment'])
        self.assertEqual(strings.share_experiment_page_title, info['pageTitle'])
        self.assertEqual(None, info['reason'])

    
    def test_my_experiments_should_show_experiments(self):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user

        exp1 = createMockExperiment(1, 0)
        exp2 = createMockExperiment(2, 0)
        
        ptm_user.myExperiments.return_value = [exp1,exp2]
        
        info = manage_experiments(request)
        
        self.assertEqual([exp1, exp2], info['experiments'])
        self.assertEqual(strings.my_experiments_page_title, info['pageTitle'])
        