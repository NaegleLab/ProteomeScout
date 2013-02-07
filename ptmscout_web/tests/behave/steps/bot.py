from ptmscout.database import user, experiment, DBSession
from mock import patch
import re
from ptmscout.database.permissions import Permission
from ptmscout.config import strings
import os
from assertions import assertContains, assertRegexMatch, synchronous_assert_called

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
        response = self.app.get('/register')
        response.form.set('username', self.username)
        response.form.set('pass1', self.password)
        response.form.set('pass2', self.password)
        response.form.set('name', self.fullname)
        response.form.set('email', self.email)
        response.form.set('institution', self.institution)
        response = response.form.submit().follow()

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
            exp = experiment.getExperimentById(exp_id, self.user)
            exp.public = 0
            self.user.permissions.append(Permission(exp, 'owner'))
        
        self.user.saveUser()
        
    def publish_experiment(self, exp_id):
        result = self.app.get("http://localhost/account/experiments")
        forms = result.forms__get()
        form_name = 'publish%d' % (exp_id)
        
        result = self.app.get(forms[form_name].action, status=200)
        result.mustcontain(strings.publish_experiment_confirm_message)
        
        forms = result.forms__get()
        result = forms['confirm'].submit()
        result.mustcontain(strings.publish_experiment_success_message)
        
    def load_datafile(self, filename, form, load_type='new', parent_experiment='', change_description=''):
        filename = os.path.join('tests','behave','data',filename)
        
        f = open(filename, 'rb')
        filecontents = f.read()
        
        
        form.set('data_file', (filename, filecontents))
        form.set('load_type', load_type)
        form.set('parent_experiment', parent_experiment)
        form.set('change_description', change_description)
        return form.submit()
    

def set_metadata_form_defaults(form):
    form.set('pmid', '')
    form.set('URL', '')
    
    form.set('description', "This is an experiment")
    
    form.set('published', "yes")
    form.set('experiment_name', "Experiment with some kind of data")
    form.set('author_contact', "author@institute.edu")
    form.set('authors', "S Guy, S O Person")
    form.set('journal', "Journal of serendipitous results")
    form.set('publication_month', "december")
    form.set('publication_year', "2011")
    form.set('volume', "236")
    form.set('page_start', "111")
    form.set('page_end', "123")
    
    form.set('ambiguous', "yes")
    form.set('notes', "none")


def upload_file(context, filename, force=False):
    context.active_user.login()
    
    context.form = context.ptmscoutapp.get('/upload',status=200).form
    context.result = context.active_user.load_datafile(filename, context.form).follow()
    context.result = context.result.form.submit()
    
    if context.result.status != '302 Found':
        if force:
            context.result.form.set('override', 'yes')
            context.result = context.result.form.submit()
        else:  
            context.result.showbrowser()
    
    assert context.result.status == '302 Found'
    context.result = context.result.follow()
    
    
    set_metadata_form_defaults(context.result.form)
    
    context.result = context.result.form.submit().follow()
    context.result = context.result.form.submit().follow()
    
    context.form = context.result.form
    context.form.set('terms_of_use', "yes")

def log_abort():
    import logging 
    logging.getLogger('ptmscout').debug("Transaction aborted")
    raise Exception()

def session_flush():
    import logging
    from ptmscout.database import DBSession
    DBSession.flush()
    logging.getLogger('ptmscout').debug("Transaction committed")

def check_experiment_loaded(context, exp_title, patch_mail):
    DBSession.flush()

    my_experiments_page = "http://localhost/account/experiments"
    result = context.result

    result.mustcontain(strings.experiment_upload_started_page_title)
    result.mustcontain(strings.experiment_upload_started_message % (my_experiments_page))
    
    result = context.ptmscoutapp.get(my_experiments_page)
    
    result.mustcontain(exp_title)
    result.mustcontain("loaded")

    m = re.search('http://localhost/experiments/([0-9]+)', str(result))
    if m == None:
        print 'http://localhost/experiments/[0-9]+ not found'
    exp_link = m.group(0)
    context.exp_id = int(m.group(1))

    exp = experiment.getExperimentById(context.exp_id)
    DBSession.refresh(exp)
    
    peps = context.table[0]['peptides']
    prots = context.table[0]['proteins']
    errors = context.table[0]['errors']
    rejected = context.table[0]['rejected']
    
    result = context.ptmscoutapp.get(exp_link + "/summary")
    
    strres = str(result).replace("\n", " ")
    assertRegexMatch('Proteins.*%s' % (prots), strres)
    assertRegexMatch('Peptides.*%s' % (peps), strres)
    assertRegexMatch('Rejected Peptides.*%s' % (rejected), strres)
    
    synchronous_assert_called(patch_mail, 5)
    argstr = str(patch_mail.call_args)

    assertContains("Peptides: %s" % (peps), argstr)
    assertContains("Proteins: %s" % (prots), argstr)
    assertContains("Errors: %s" % (errors), argstr)
    assertRegexMatch('http://localhost/experiments/[0-9]+/errors', argstr)
    
    context.exp_link = exp_link
    
    browse_page = context.exp_link + "/browse"
    result = context.ptmscoutapp.get(browse_page)
    context.result = result
