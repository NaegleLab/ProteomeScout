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
    

@given(u'a user submits a dataset with an accession that looks like a GenPept accession, but is not a valid accession number')
def impl(context):
    assert False

@given(u'a user has loaded a dataset and the modification type does not match the amino acid for that species')
def impl(context):
    assert False

@given(u'a user submits a correctly formatted dataset of phosphorylation data')
def impl(context):
    assert False

@given(u'a user has loaded a dataset in which modifications have varying specificities of naming')
def impl(context):
    assert False

@given(u'a user submits a dataset in which a single peptide has more than one flavor of modification and the mod_type entry does not match the amino acids')
def impl(context):
    assert False

@given(u'a user submits a dataset in which an isoform specific record is included')
def impl(context):
    assert False

@given(u'a user submits a dataset that has an incorrect peptide to protein match')
def impl(context):
    assert False

@given(u'a user submits a dataset and a single peptide has more than one flavor of modification and mod_type assignments match the sequential order of modifcations')
def impl(context):
    assert False


    
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
