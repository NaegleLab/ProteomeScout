from behave import *
import re
from tests.behave.steps import bot
import os
from mock import patch
from ptmscout.database import DBSession, experiment

def log_abort():
    import logging 
    logging.getLogger('ptmscout').debug("Transaction aborted")

def session_flush():
    import logging 
    DBSession.flush()
    logging.getLogger('ptmscout').debug("Transaction committed")

@given(u'a user uploads a dataset with errors')
@patch('transaction.abort')
@patch('transaction.commit')
@patch('ptmscout.utils.mail.celery_send_mail')
def user_upload_with_errors(context, patch_mail, patch_commit, patch_abort):
    
    patch_commit.side_effect = session_flush
    patch_abort.side_effect = log_abort
    
    context.active_user.login()
    bot.upload_file(context, 'datasetLoad_badAminoAcid.txt', True)
    context.form.submit()

    
    result = context.ptmscoutapp.get("http://localhost/account/experiments")
    
    m = re.search('http://localhost/experiments/([0-9]+)', str(result))
    if m == None:
        print 'http://localhost/experiments/([0-9]+) not found'
    context.exp_link = m.group(0)
    context.exp_id = int(m.group(1))

    DBSession.flush()
    exp = experiment.getExperimentById(context.exp_id, secure=False)
    DBSession.refresh(exp)
    
    context.active_user.publish_experiment(context.exp_id)
    
@when(u'the user opens the error log')
def user_opens_the_error_log(context):
    context.result = context.ptmscoutapp.get(context.exp_link + "/errors")

@then(u'the user should see a list of the proteins and peptides that were rejected during data loading')
def show_proteins_and_rejected_peptides(context):

    
    context.result.mustcontain('Q6XBG2')
    context.result.mustcontain('P60711')
    
    context.result.mustcontain('YLRFSDIkKNINSGA')
    context.result.mustcontain('TLkYPIEhGIVTNWD')
    context.result.mustcontain('TLKYPIEhgIVTNWD')
    
    context.result.mustcontain('6')
    context.result.mustcontain('8')
    context.result.mustcontain('9')
    
    context.result.mustcontain("Warning: Specified modification 'lysine methyl ester' does not match residue 'K' for specified species")
    context.result.mustcontain("Warning: Specified modification 'methylhistidine' does not match residue 'K' for any known species")
    context.result.mustcontain("Warning: Specified modification 'methylation' does not match residue 'G' for any known species")


@then(u'the user should be given a link to download the subset of the dataset that did not load with an additional error explanation column')
def download_annotated_dataset(context):
    download_link = 'http://localhost/experiments/%d/download?errors=True' % (context.exp_id)
    context.result.mustcontain(download_link)
    
    result = context.ptmscoutapp.get(download_link)
    
    expected_result = open(os.path.join('tests','behave','data','datasetLoad_errorAnnotations.txt')).read()
    result.mustcontain(expected_result)

@then(u'the user should be given a link to append data to the dataset')
def append_to_dataset(context):
    append_link = 'http://localhost/upload?load_type=append&parent_experiment=%d' % (context.exp_id)
    context.result.mustcontain(append_link)
    
    result = context.ptmscoutapp.get(append_link)
    
    result.mustcontain('<option value="%d" selected>Experiment with some kind of data</option>' % (context.exp_id))
    result.mustcontain('name="load_type" value="append" checked />')
    
    
