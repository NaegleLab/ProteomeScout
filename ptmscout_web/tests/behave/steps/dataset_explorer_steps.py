from behave import *
import bot
from ptmscout.config import settings, strings
from mock import patch
from tests.behave.steps.assertions import assertContains, assertEqual, assertIn,\
    assertAlmostEqual
import os
import base64
import json


@given(u'a user uploads annotations to an experiment')
def user_upload_annotations(context):
    context.active_user.login()
    
    result = context.ptmscoutapp.get('/experiments/28/subsets')
    annote_file = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'experiment.28.test.annotations.tsv')
    bot.set_file_form_contents('annotationfile', annote_file, result.forms[1])
    
    
    result = result.forms[1].submit(status=302)
    result = result.follow(status=200)
    
    result = result.form.submit()
    if result.status != '302 Found':
        result.form.set('override', 'yes')
        result = result.form.submit()
    
    result = result.follow()
    
    
    result.form.set('terms_of_use', "yes")
    context.result = result


@then(u'the user should be sent an email with a link to the dataset explorer with')
@patch('transaction.abort')
@patch('transaction.commit')
@patch('ptmscout.utils.mail.celery_send_mail')
def user_receives_email(context, patch_mail, patch_commit, patch_abort):
    patch_commit.side_effect = bot.session_flush
    patch_abort.side_effect = bot.log_abort
    
    context.result.form.submit()    
    
    argstr = str(patch_mail.call_args)
    
    assertContains(context.active_user.user.email, argstr)
    assertContains(strings.annotation_upload_finished_subject, argstr)
    assertContains("Annotations Loaded: %s" % (context.table[0]['annotations']), argstr)
    assertContains("Errors Encountered: %s" % (context.table[0]['errors']), argstr)
    
    
@then(u'the user should be able to see their annotations in the dataset explorer')
def user_views_dataset_explorer(context):
    result = context.ptmscoutapp.get("/experiments/28/subsets")
    
    d = result.pyquery
    
    json_base64 = d('#field-data').text()
    response = json.loads( base64.b64decode(json_base64) )
    
    assertIn( 'Fun', response['metadata_fields'] )
    assertEqual(['fun', 'not fun'], response['metadata_fields']['Fun'])
    
    assertIn( 'Random', response['clustering_labels'] )
    assertEqual( ['1','2','3','4','5','6'], response['clustering_labels']['Random'] )
    
    assertIn( 'Interestingness', response['quantitative_fields'] )
    
@then(u'the user should be able to view the clusterings')
def user_should_be_able_to_query_clusters(context):
    for i in xrange(1, 7):
        fq = [ 
              ['nop', ['cluster', 'Random', 'eq', str(i)]], 
             ]
        result = context.active_user.query_dataset_explorer(28, fq)
        assert 9 <= result['foreground']['peptides']

@then(u'the user should be able to filter by numerical data')
def filter_ms_with_greater_than_annotated_value(context):
    fq = [
          ['nop', ['quantitative', 'Interestingness', 'gt', '0.5']]
         ]
    
    result = context.active_user.query_dataset_explorer( 28, fq )
    assertEqual( 26, result['foreground']['peptides'] )
    
    context.result = result
    
@then(u'the user should be able to filter by nominative data')
def filter_nominative_fun(context):
    fq = [ 
          ['nop', ['metadata', 'Fun', 'eq', 'fun']], 
         ]
    result = context.active_user.query_dataset_explorer(28, fq)
    assertEqual(6, result['foreground']['peptides'])

    
@then(u'the user should be able to view nominative features in the enrichment data')
def view_nominative_enrichment(context):
    enrichment_labels = set( [ row[0] for row in context.result['enrichment'] ] )
    enrichment_values = set( [ row[1] for row in context.result['enrichment'] ] )
    
    fun_p_values = {}
    for row in context.result['enrichment']:
        if row[0] == "Annotation: Fun":
            fun_p_values[row[1]] = row[2]
    
    assertIn('Annotation: Fun', enrichment_labels)
    assertIn('fun', enrichment_values)
    assertIn('not fun', enrichment_values)
    
    assertAlmostEqual(0.23759270964480667, fun_p_values['not fun'])
    assertAlmostEqual(0.3949453933540862, fun_p_values['fun'])



@given(u'a user uploads a file containing non-mass spec experimental data')
def upload_dataset_for_annotation(context):
    pass
    
    