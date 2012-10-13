from behave import *
from ptmscout.config import strings
import time
from assertions import assertContains
from mock import patch

def set_form_defaults(form):
    form.set('pmid', '')
    form.set('URL', '')
    
    form.set('published', "yes")
    form.set('experiment_name', "Experiment with correct data")
    form.set('author_contact', "author@institute.edu")
    form.set('authors', "S Guy, S O Person")
    form.set('journal', "Journal of serendipitous results")
    form.set('publication_month', "december")
    form.set('publication_year', "2011")
    form.set('volume', "236")
    form.set('page_start', "111")
    form.set('page_end', "123")
    
    form.set('ambiguous', "no")
    form.set('notes', "none")
    form.set('terms_of_use', "yes")


@given(u'a user submits a correctly formatted dataset of phosphorylation data')
def submit_correct_dataset(context):
    context.active_user.login()
    result = context.ptmscoutapp.get('/upload', status=200)
    context.form = result.form
    
    context.form.set('load_type', "new")
    
    filename = "tests/behave/data/datasetLoad_correctDataset.txt"
    f = open(filename, 'rb')
    filecontents = f.read()
    
    context.form.set('data_file', (filename, filecontents))
    context.form.set('description', "This is a correct dataset")
    
    set_form_defaults(context.form)

@given(u'a user submits a dataset that has an incorrect peptide to protein match')
def submit_incorrect_peptide_dataset(context):
    context.active_user.login()
    result = context.ptmscoutapp.get('/upload', status=200)
    context.form = result.form
    
    context.form.set('load_type', "new")
    filename = "tests/behave/data/datasetLoad_pepMismatch.txt"
    f = open(filename, 'rb')
    filecontents = f.read()
    
    context.form.set('data_file', (filename, filecontents))
    context.form.set('description', "This is a dataset with a bad peptide")
    
    set_form_defaults(context.form)
    

@given(u'a user submits a dataset that has an accession that looks like a GenPept accession, but is not currently there')
def submit_incorrect_genpept_dataset(context):
    context.active_user.login()
    result = context.ptmscoutapp.get('/upload', status=200)
    context.form = result.form
    
    context.form.set('load_type', "new")
    filename = "tests/behave/data/datasetLoad_badAcc.txt"
    f = open(filename, 'rb')
    filecontents = f.read()
    
    context.form.set('data_file', (filename, filecontents))
    context.form.set('description', "This is a dataset with a bad accession")
    
    set_form_defaults(context.form)
    

@given(u'a user submits a dataset file with more than one peptide column')
def submit_dataset_with_multiple_peptide_columns(context):
    context.active_user.login()
    result = context.ptmscoutapp.get('/upload', status=200)
    context.form = result.form
    
    context.form.set('load_type', "new")
    filename = "tests/behave/data/datasetLoad_moreThanOnePep.txt"
    f = open(filename, 'rb')
    filecontents = f.read()
    
    context.form.set('data_file', (filename, filecontents))
    context.form.set('description', "This is a dataset with multiple peptide columns")
    
    set_form_defaults(context.form)
    

@given(u'a user submits a dataset with no accession column')
def submit_dataset_without_accession_column(context):
    context.active_user.login()
    result = context.ptmscoutapp.get('/upload', status=200)
    context.form = result.form
    
    context.form.set('load_type', "new")
    filename = "tests/behave/data/datasetLoad_noAcc.txt"
    f = open(filename, 'rb')
    filecontents = f.read()
    
    context.form.set('data_file', (filename, filecontents))
    context.form.set('description', "This is a dataset without an accession column")
    
    set_form_defaults(context.form)
    

@then(u'the user should see an error that says "{error_text}"')
def error_text_display(context, error_text):
    context.result = context.form.submit()
    context.result.mustcontain(error_text)
    
    
@then(u'the user should see the text "{text}"')
def info_text_display(context, text):
    context.result.mustcontain(text)
    
@then(u'the user should be sent an email with a link to the experiment which contains')
@patch('ptmscout.utils.mail.send_automail_message')
def experiment_uploaded_check_email(context, patch_mail):
    result = context.form.submit().follow()

    #confirm the upload
    result = result.form.submit()
    
    result.mustcontain(strings.experiment_upload_started_page_title)
    result.mustcontain(strings.experiment_upload_started_message[0:-len("<a href=\"%s\">this page</a>")])
    
    asynchronous_assert_called(patch_mail)
    argstr = str(patch_mail.call_args)
    
    peps = context.table['peptides']
    prots = context.table['proteins']
    errors = context.table['errors']
    
    assertContains("Peptides: %s" % (peps), argstr)
    assertContains("Proteins: %s" % (prots), argstr)
    assertContains("Errors: %s" % (errors), argstr)
    

def asynchronous_assert_called(mock, limit=1):
    slept_time = 0.0
    while not mock.called and slept_time < limit:
        time.sleep(0.1)
        slept_time += 0.1    
    
    assert mock.called, "Asynchronous call to %s timed out" % (str(mock))
