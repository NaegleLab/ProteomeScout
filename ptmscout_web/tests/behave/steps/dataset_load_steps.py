from behave import *
from ptmscout.config import strings
import time
from assertions import assertContains, assertRegexMatch
from mock import patch
import re
from tests.behave.steps import bot


@given(u'a user submits a dataset with an accession that looks like a GenPept accession, but is not a valid accession number')
def submit_bad_accession(context):
    bot.upload_file(context, 'datasetLoad_badAcc.txt')

@given(u'a user has loaded a dataset and the modification type does not match the amino acid for that species')
def submit_bad_modification_type_for_amino_acid(context):
    bot.upload_file(context, 'datasetLoad_badAminoAcid.txt', True)

@given(u'a user submits a correctly formatted dataset of phosphorylation data')
def submit_correct_dataset(context):
    bot.upload_file(context, 'datasetLoad_correctDataset.txt')

@given(u'a user has loaded a dataset in which modifications have varying specificities of naming')
def submit_dataset_with_degenerate_names(context):
    bot.upload_file(context, 'datasetLoad_methylation.txt')

@given(u'a user submits a dataset in which a single peptide has more than one flavor of modification and the mod_type entry does not match the amino acids')
def submit_mod_type_mismatch_dataset_multiple_mods(context):
    bot.upload_file(context, 'datasetLoad_moreThanOneModWOAssignment.txt', True)

@given(u'a user submits a dataset in which an isoform specific record is included')
def submit_isoform_dataset(context):
    bot.upload_file(context, 'datasetLoad_isoform.txt')

@given(u'a user submits a dataset that has an incorrect peptide to protein match')
def submit_incorrect_peptide_dataset(context):
    bot.upload_file(context, 'datasetLoad_pepMismatch.txt')

@given(u'a user submits a dataset and a single peptide has more than one flavor of modification and mod_type assignments match the sequential order of modifcations')
def submit_multiple_mods_correct_dataset(context):
    bot.upload_file(context, 'datasetLoad_moreThanOneMod.txt')

def log_abort():
    import logging 
    logging.getLogger('ptmscout').debug("Transaction aborted")

def session_flush():
    import logging 
    from ptmscout.database import DBSession
    DBSession.flush()
    logging.getLogger('ptmscout').debug("Transaction committed")
    
@then(u'the user should be sent an email with a link to the experiment which contains')
@patch('transaction.abort')
@patch('transaction.commit')
@patch('ptmscout.utils.mail.celery_send_mail')
def experiment_uploaded_check_email(context, patch_mail, patch_commit, patch_abort):
    patch_commit.side_effect = session_flush
    patch_abort.side_effect = log_abort
    
    result = context.form.submit()
    
    my_experiments_page = "http://localhost/account/experiments"
    
    result.mustcontain(strings.experiment_upload_started_page_title)
    result.mustcontain(strings.experiment_upload_started_message % (my_experiments_page))
    
    result = context.ptmscoutapp.get(my_experiments_page)
    
    result.mustcontain("Experiment with some kind of data")
    result.mustcontain("loaded")
    
    m = re.search('http://localhost/experiments/([0-9]+)', str(result))
    if m == None:
        print 'http://localhost/experiments/[0-9]+ not found'
    exp_link = m.group(0)
    context.exp_id = int(m.group(1))
    
    peps = context.table[0]['peptides']
    prots = context.table[0]['proteins']
    errors = context.table[0]['errors']
    
    result = context.ptmscoutapp.get(exp_link + "/summary")
    
    strres = str(result).replace("\n", " ")
    assertRegexMatch('Proteins.*%s' % (prots), strres)
    assertRegexMatch('Peptides.*%s' % (peps), strres)
    assertRegexMatch('Rejected Peptides.*%s' % (errors), strres)
    
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

@then(u'the experiment browser should contain the correct peptides')
def check_displayed_peptides(context):
    result = context.result
    
    strres = str(result).replace("\n", " ")
    #       1390       1400       1410
    # AGLGIRQGGK APVTPRGRGR RGRPPSRTTG
    #   LGIRQGGk APVTPRG
    #         GK APVTPrGRGR RGR
    #            APVTPRGrGR RGRPP
    #              VTPRGRGr RGRPPSR
    #               TPRGRGR rGRPPSRT
    #                 RGRGR RGrPPSRTTG
    
    assertRegexMatch("%s.*%s" % ("K1390", "LGIRQGGkAPVTPRG"), strres)
    assertRegexMatch("%s.*%s" % ("R1396", "GKAPVTPrGRGRRGR"), strres)
    assertRegexMatch("%s.*%s" % ("R1398", "APVTPRGrGRRGRPP"), strres)
    assertRegexMatch("%s.*%s" % ("R1400", "VTPRGRGrRGRPPSR"), strres)
    assertRegexMatch("%s.*%s" % ("R1401", "TPRGRGRrGRPPSRT"), strres)
    assertRegexMatch("%s.*%s" % ("R1403", "RGRGRRGrPPSRTTG"), strres)
    
    assertRegexMatch("%s.*%s" % ("K317", "FSFKLLKkECPIPNV"), strres)
    assertRegexMatch("%s.*%s" % ("H73", "TLKYPIEhGIVTNWD"), strres)
    assertRegexMatch("%s.*%s" % ("R144", "ENFVDKLrESLMSVA"), strres)
    assertRegexMatch("%s.*%s" % ("H73", "TLKYPIEhGIVTNWD"), strres)
    assertRegexMatch("%s.*%s" % ("R257", "GTGPQRPrSWAAADS"), strres)
    assertRegexMatch("%s.*%s" % ("K212", "TLNEDSYkDSTLIMQ"), strres)
    

def synchronous_assert_called(mock, limit=1):
    slept_time = 0.0
    while not mock.called and slept_time < limit:
        time.sleep(0.1)
        slept_time += 0.1    
    
    assert mock.called, "Synchronous call to %s timed out" % (str(mock))
