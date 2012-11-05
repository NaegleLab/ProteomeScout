from behave import *
from ptmscout.config import strings
import time
from assertions import assertContains
from mock import patch

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
    
    form.set('ambiguous', "no")
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

@given(u'a user submits a dataset with an accession that looks like a GenPept accession, but is not a valid accession number')
def submit_bad_accession(context):
    upload_file(context, 'datasetLoad_badAcc.txt')

@given(u'a user has loaded a dataset and the modification type does not match the amino acid for that species')
def submit_bad_modification_type_for_amino_acid(context):
    upload_file(context, 'datasetLoad_badAminoAcid.txt', True)

@given(u'a user submits a correctly formatted dataset of phosphorylation data')
def submit_correct_dataset(context):
    upload_file(context, 'datasetLoad_correctDataset.txt')

@given(u'a user has loaded a dataset in which modifications have varying specificities of naming')
def submit_dataset_with_degenerate_names(context):
    upload_file(context, 'datasetLoad_methylation.txt')

@given(u'a user submits a dataset in which a single peptide has more than one flavor of modification and the mod_type entry does not match the amino acids')
def submit_mod_type_mismatch_dataset_multiple_mods(context):
    upload_file(context, 'datasetLoad_moreThanOneModWOAssignment.txt', True)

@given(u'a user submits a dataset in which an isoform specific record is included')
def submit_isoform_dataset(context):
    upload_file(context, 'datasetLoad_isoform.txt')

@given(u'a user submits a dataset that has an incorrect peptide to protein match')
def submit_incorrect_peptide_dataset(context):
    upload_file(context, 'datasetLoad_pepMismatch.txt')

@given(u'a user submits a dataset and a single peptide has more than one flavor of modification and mod_type assignments match the sequential order of modifcations')
def submit_multiple_mods_correct_dataset(context):
    upload_file(context, 'datasetLoad_moreThanOneMod.txt')


    
@then(u'the user should be sent an email with a link to the experiment which contains')
@patch('ptmscout.utils.mail.send_automail_message')
def experiment_uploaded_check_email(context, patch_mail):
    result = context.form.submit()
    
    result.mustcontain(strings.experiment_upload_started_page_title)
    result.mustcontain(strings.experiment_upload_started_message[0:-len("<a href=\"%s\">this page</a>")])
    
    synchronous_assert_called(patch_mail)
    argstr = str(patch_mail.call_args)
    
    peps = context.table[0]['peptides']
    prots = context.table[0]['proteins']
    errors = context.table[0]['errors']
    
    assertContains("Peptides: %s" % (peps), argstr)
    assertContains("Proteins: %s" % (prots), argstr)
    assertContains("Errors: %s" % (errors), argstr)
    

def synchronous_assert_called(mock, limit=1):
    slept_time = 0.0
    while not mock.called and slept_time < limit:
        time.sleep(0.1)
        slept_time += 0.1    
    
    assert mock.called, "Synchronous call to %s timed out" % (str(mock))
