from behave import *
from ptmscout.config import strings
import time
from assertions import assertContains, assertRegexMatch
from mock import patch
import re
from tests.behave.steps import bot
import helpers, assertions

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

@given(u'a user submits a dataset in which no protein accession is found')
def submit_all_error_dataset(context):
    bot.upload_file(context, 'datasetLoad_allFail.txt', True)

@given(u'a user submits a dataset with 1000 measured peptides')
def submit_big_dataset(context):
    bot.upload_file(context, 'datasetLoad_big_dataset.txt')

@given(u'a user submits a dataset in which viral proteins are included')
def submit_viral_dataset(context):
    bot.upload_file(context, 'datasetLoad_virus.txt')

@given(u'a user submits a dataset in which the same peptide')
def submit_dataset_with_repeat_peptides_diff_mod(context):
    bot.upload_file(context, 'datasetLoad_samePepDiffMod.txt')

   
@then(u'the user should be sent an email with a link to the experiment which contains')
@patch('transaction.abort')
@patch('transaction.commit')
@patch('ptmscout.utils.mail.celery_send_mail')
def experiment_uploaded_check_email(context, patch_mail, patch_commit, patch_abort):
    patch_commit.side_effect = bot.session_flush
    patch_abort.side_effect = bot.log_abort

    context.result = context.form.submit()

    bot.check_experiment_loaded(context, 'Experiment with some kind of data', patch_mail)


@then(u'the experiment browser should contain the correct peptides')
def check_displayed_peptides(context):
    p = context.result.pyquery
    
    table_tag = p('table')
    pytable = helpers.parse_table(p, table_tag)
    
    exp_table = \
    [['Protein', 'Gene', 'Species', 'Sequence Length', '#Reported Sources', '#Modified Residues', 'Modified Amino Acids', 'Modification Types'],
     ['AP2-associated protein kinase 1', 'AAK1', 'homo sapiens', '961', '1', '1', 'K', 'Trimethylation'],
     ['4-aminobutyrate aminotransferase, mitochondrial', 'Abat', 'rattus norvegicus', '500', '1', '1', 'R', 'Asymmetric dimethylarginine'],
     ['Actin, cytoplasmic 1; AltName: Full=Beta-actin; Contains: Actin, cytopla', 'Actb', 'rattus norvegicus', '375', '2', '1', 'H', 'Phosphohistidine'],
     ['Actin, cytoplasmic 1; AltName: Full=Beta-actin; Contains: Actin, cytopla', 'ACTB', 'homo sapiens', '375', '2', '1', 'H', 'Methylhistidine'],
     ['Protein FAM228A', 'Fam228a', 'mus musculus', '299', '1', '1', 'R', 'Dimethylation'],
     ['Tumor suppressor p53-binding protein 1; Short=p53-binding protein 1; Short=p53BP1; Sho', 'TP53BP1', 'homo sapiens', '1972', '2', '6', 'K,R', 'N6-methyllysine, Omega-N-methylarginine'],
     ['14-3-3 protein theta; AltName: Full=14-3-3 protein tau; AltName: Full=14-3-3 protein T', 'YWHAQ', 'homo sapiens', '245', '2', '1', 'K', 'N6-methyllysine']]

    assertions.assertEqual(8, len(pytable))
    for i in xrange(0, 8):
        assertions.assertIn(exp_table[i], pytable)


