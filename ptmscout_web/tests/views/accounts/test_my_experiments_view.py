from mock import patch, call
from ptmscout.config import strings
from ptmscout.database.user import NoSuchUser
from ptmscout.views.accounts.my_experiments_view import confirm_invite_user, \
    privatize_experiment, publish_experiment, manage_experiment_permissions, \
    manage_experiments, get_sessions, remove_user_permission
from pyramid.httpexceptions import HTTPFound
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockUser, createMockExperiment, \
    createMockPermission, createMockSession, createMockJob
from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.database import upload

class MyExperimentsViewIntegrationTests(IntegrationTestCase):
    
    def test_integration_with_datasets(self):
        from ptmscout.database import experiment
        
        exp = experiment.getExperimentById(26, None, False)
        exp.type='dataset'
        exp.job.status = 'running'
        exp.job.stage = 'proteins'
        exp.job.progress = 1000
        exp.job.max_progress = 10000
        exp.job.user_id = self.bot.user.id
        exp.grantPermission(self.bot.user, 'owner')
        exp.saveExperiment()

        exp2 = experiment.getExperimentById(25, None, False)
        exp2.type='dataset'
        exp2.job.status = 'running'
        exp2.job.stage = 'GO terms'
        exp2.job.progress = 0
        exp2.job.max_progress = 0
        exp2.job.user_id = self.bot.user.id
        exp2.grantPermission(self.bot.user, 'owner')
        exp2.saveExperiment()


        exp3 = experiment.getExperimentById(28, None, False)
        exp3.type='dataset'
        exp3.job_id = None
        exp3.grantPermission(self.bot.user, 'owner')
        exp3.saveExperiment()

        exp4 = experiment.getExperimentById(1, None, False)
        exp4.type='dataset'
        exp4.job.status = 'error'
        exp4.job.user_id = self.bot.user.id
        exp4.grantPermission(self.bot.user, 'owner')
        exp4.saveExperiment()

        session = upload.Session()
        session.experiment_id = exp3.id
        session.user_id = self.bot.user.id
        session.data_file = ''
        session.change_name = ''
        session.change_description = ''
        session.load_type = 'new'
        session.resource_type = 'dataset'
        session.stage = 'confirm'
        session.save()

        session2 = upload.Session()
        session2.experiment_id = exp4.id
        session2.user_id = self.bot.user.id
        session2.data_file = ''
        session2.change_name = ''
        session2.change_description = ''
        session2.load_type = 'new'
        session2.resource_type = 'dataset'
        session2.stage = 'confirm'
        session2.save()
        
        self.bot.login()
        
        result = self.ptmscoutapp.get('/account/experiments')
        
        result.mustcontain(exp.name)
        result.mustcontain('Status')
        result.mustcontain('processing: proteins')
        result.mustcontain('1000 / 10000')
        result.mustcontain('processing: GO terms')
        result.mustcontain('error')
        result.mustcontain('retry')
        result.mustcontain('continue upload')
        
        result.mustcontain('<a href="%s/dataset/upload/%d">continue upload</a>' % ("http://localhost", session.id))
        result.mustcontain('<a href="%s/dataset/upload/%d/retry">retry</a>' % ("http://localhost", session2.id))
    
    def test_integration(self):
        from ptmscout.database import experiment
        
        exp = experiment.getExperimentById(26, None, False)
        exp.job.status = 'running'
        exp.job.stage = 'proteins'
        exp.job.progress = 1000
        exp.job.max_progress = 10000
        exp.job.user_id = self.bot.user.id
        exp.grantPermission(self.bot.user, 'owner')
        exp.saveExperiment()

        exp2 = experiment.getExperimentById(25, None, False)
        exp2.job.status = 'running'
        exp2.job.stage = 'GO terms'
        exp2.job.progress = 0
        exp2.job.max_progress = 0
        exp2.job.user_id = self.bot.user.id
        exp2.grantPermission(self.bot.user, 'owner')
        exp2.saveExperiment()


        exp3 = experiment.getExperimentById(28, None, False)
        exp3.job_id = None
        exp3.grantPermission(self.bot.user, 'owner')
        exp3.saveExperiment()

        exp4 = experiment.getExperimentById(1, None, False)
        exp4.job.status = 'error'
        exp4.job.user_id = self.bot.user.id
        exp4.grantPermission(self.bot.user, 'owner')
        exp4.saveExperiment()

        session = upload.Session()
        session.experiment_id = exp3.id
        session.user_id = self.bot.user.id
        session.data_file = ''
        session.change_name = ''
        session.change_description = ''
        session.load_type = 'new'
        session.resource_type = 'experiment'
        session.stage = 'confirm'
        session.save()

        session2 = upload.Session()
        session2.experiment_id = exp4.id
        session2.user_id = self.bot.user.id
        session2.data_file = ''
        session2.change_name = ''
        session2.change_description = ''
        session2.load_type = 'new'
        session2.resource_type = 'experiment'
        session2.stage = 'confirm'
        session2.save()
        
        self.bot.login()
        
        result = self.ptmscoutapp.get('/account/experiments')
        
        result.mustcontain(exp.name)
        result.mustcontain('Status')
        result.mustcontain('processing: proteins')
        result.mustcontain('1000 / 10000')
        result.mustcontain('processing: GO terms')
        result.mustcontain('error')
        result.mustcontain('retry')
        result.mustcontain('continue upload')
        
        result.mustcontain('<a href="%s/upload/%d">continue upload</a>' % ("http://localhost", session.id))
        result.mustcontain('<a href="%s/upload/%d/retry">retry</a>' % ("http://localhost", session2.id))
        

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
   
    def test_remove_user_permission_should_delete_correct_user(self):
        u1 = createMockUser()
        u2 = createMockUser()
        u3 = createMockUser()

        request = DummyRequest()
        request.POST['submitted'] = "1"
        request.POST['mode'] = 'add'
        request.POST['uid'] = str(u2.id)

        users = [u1,u2,u3]
        exp = createMockExperiment()

        remove_user_permission(request, users, exp)

        exp.revokePermission.assert_called_once_with(u2)
        exp.saveExperiment.assert_called_once_with()
        self.assertEqual([u1,u3], users)

    @patch('ptmscout.views.accounts.my_experiments_view.remove_user_permission')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_manage_experiment_permissions_should_call_remove_on_remove_user(self, patch_getExperiment, patch_remove):
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
        request.POST['mode'] = 'remove'
        request.POST['email'] = "100"

        patch_remove.return_value = None

        result = manage_experiment_permissions(request)

        self.assertEqual(None, result['reason'])
        self.assertEqual(exp1, result['experiment'])
        self.assertEqual([user1], result['users'])
        self.assertEqual(strings.share_experiment_page_title, result['pageTitle'])

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
        request.POST['mode'] = 'add'
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
        request.POST['mode'] = 'add'
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
    
    
    @patch('ptmscout.database.upload.getMostRecentSession')
    def test_get_sessions_should_get_most_recent_sessions_for_experiments(self, patch_getSession):
        user = createMockUser()
        exp = createMockExperiment(10)
        session = createMockSession(user)
        
        patch_getSession.return_value = session
        
        smap = get_sessions([exp])
        
        patch_getSession.assert_called_once_with(10)
        
        self.assertEqual({exp.id: session.id}, smap)

    @patch('ptmscout.views.accounts.my_experiments_view.get_sessions')
    def test_my_experiments_should_show_experiments(self, patch_getSessions):
        request = DummyRequest()
        ptm_user = createMockUser("username", "email", "password", 1)
        request.user = ptm_user
        request.user.jobs = [ createMockJob(), createMockJob(), createMockJob() ]

        exp1 = createMockExperiment(1, 0)
        exp2 = createMockExperiment(2, 0)
        exp3 = createMockExperiment(3, 0)
        exp4 = createMockExperiment(4, 0)
        exp1.status=  'finished'
        exp2.status = 'running'
        exp3.status = 'configuration'
        exp4.status = 'error'
        
        patch_getSessions.return_value = {"some map":"of session ids"}
        
        ptm_user.myExperiments.return_value = [exp1, exp2, exp3, exp4]
        ptm_user.myDatasets.return_value = []

        info = manage_experiments(request)
        
        self.assertTrue( call([exp3, exp4]) in patch_getSessions.call_args_list)

        self.assertEqual([exp1, exp2, exp4], info['experiments']['available'])
        self.assertEqual([exp3], info['experiments']['in_process'])
        self.assertEqual(patch_getSessions.return_value, info['experiments']['sessions'])
        
        self.assertEqual([], info['datasets']['available'])
        self.assertEqual([], info['datasets']['in_process'])
        self.assertEqual(patch_getSessions.return_value, info['datasets']['sessions'])
        
        self.assertEqual(strings.my_experiments_page_title, info['pageTitle'])
        
