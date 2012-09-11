from ptmscout.database import user, experiment
from mock import patch
import re
from ptmscout.database.permissions import Permission
from ptmscout import strings

class Bot(object):
    def __init__(self, app):
        from ptmscout.utils.crypto import randomString
        self.app = app
        self.ext = randomString(5)
        self.username = "user_%s" % (self.ext)
        self.password = "pass%s" % (self.ext)
        self.fullname = "User %s" % (self.ext)
        self.email = "user%s@institute.edu" % (self.ext)
        self.institution = "The Institute of Turtles"

    def register_and_login(self):
        self.register()
        self.activate()
        self.login()

    # this needs to change to a direct request for registration
    @patch('ptmscout.utils.mail.send_automail_message')
    def register(self, patch_mail):
        response = self.app.post('/process_registration', {'username':self.username, 'pass1':self.password, 'pass2':self.password, 'name':self.fullname, 'email':self.email, 'institution':self.institution}, status=200)
        response.mustcontain(strings.user_registration_success_header)

        patch_mail.assert_called_once()
        self.activation_token = re.compile("&token=([a-z0-9]+)\">this link").search(patch_mail.call_args.__repr__()).group(1)
        self.activation_link = re.compile("a href=\"(.*)\">this link").search(patch_mail.call_args.__repr__()).group(1)

        self.user = user.getUserByUsername(self.username)
        
    def activate(self):
        response = self.app.get('/activate_account', {'username':self.username, 'token':self.activation_token}, status=200)
        response.mustcontain(strings.account_activation_success_header)
        
    def login(self):
        response = self.app.post('/process_login', {'username':self.username, 'password':self.password}, status=200)
        response.mustcontain(strings.login_page_success_message)
        
        self.user = user.getUserByUsername(self.username)
        
    def logout(self):
        response = self.app.get('/logout', status=200)
        response.mustcontain(strings.logout_page_message)
        
        self.user = None
        
    def acquire_experiments(self, exp_ids):
        for exp_id in exp_ids:    
            exp = experiment.getExperimentById(exp_id)
            exp.public = 0
            self.user.permissions.append(Permission(exp, 'owner'))
        
        self.user.saveUser()
        
        for exp_id in exp_ids:
            self.user.makeOwner(exp_id)
        
        