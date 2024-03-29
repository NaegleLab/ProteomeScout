from behave import *
import bot
from ptmscout.config import settings, strings
from mock import patch
from tests.behave.steps.assertions import assertContains, assertEqual, assertIn,\
    assertAlmostEqual, assertRegexMatch, synchronous_assert_called
import os
import base64
import json
import re
from ptmscout.database import experiment
import csv
import zipfile
import shutil
from tests.behave.steps import assertions


@given(u'a user uploads annotations to an experiment')
def user_upload_annotations(context):
    context.active_user.login()
    
    result = context.ptmscoutapp.get('/experiments/28/subsets')
    annote_file = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'experiment.28.test.annotations.tsv')
    result.forms[1].set('annotationname', "Annotation Set 1")
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

    annotation_form = result.forms['select-annotation-form']
    context.annotation_set_id = int(annotation_form['annotation_set'].options[1][0])
    annotation_form['annotation_set'] = str(context.annotation_set_id)

    result = annotation_form.submit()
    
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
        result = context.active_user.query_dataset_explorer(28, context.annotation_set_id, fq)
        assert 9 <= result['foreground']['peptides']

@then(u'the user should be able to filter by numerical data')
def filter_ms_with_greater_than_annotated_value(context):
    fq = [
          ['nop', ['quantitative', 'Interestingness', 'gt', '0.5']]
         ]
    
    result = context.active_user.query_dataset_explorer( 28, context.annotation_set_id, fq )
    assertEqual( 26, result['foreground']['peptides'] )
    
    context.result = result
    
@then(u'the user should be able to filter by nominative data')
def filter_nominative_fun(context):
    fq = [ 
          ['nop', ['metadata', 'Fun', 'eq', 'fun']], 
         ]
    result = context.active_user.query_dataset_explorer(28, context.annotation_set_id, fq)
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

@given(u'the user saves a subset selection from that experiment')
@patch('transaction.abort')
@patch('transaction.commit')
@patch('ptmscout.utils.mail.celery_send_mail')
def save_subset_query_using_annotation_set(context, patch_mail, patch_commit, patch_abort):
    patch_commit.side_effect = bot.session_flush
    patch_abort.side_effect = bot.log_abort
    
    context.result.form.submit()
    
    argstr = str(patch_mail.call_args)
    
    assertContains(context.active_user.user.email, argstr)
    assertContains(strings.annotation_upload_finished_subject, argstr)

    result = context.ptmscoutapp.get("/experiments/28/subsets")

    annotation_form = result.forms['select-annotation-form']
    context.annotation_set_id = int(annotation_form['annotation_set'].options[1][0])

    fq = [
          ['nop', ['quantitative', 'Interestingness', 'gt', '0.5']]
         ]
    
    result = context.active_user.save_subset_selection( 28, context.annotation_set_id, fq )
    assertEqual( 'success', result['status'])
    context.result = result

@when(u'the user creates a sharing token for the experiment subset')
def create_share_token_for_subset(context):
    context.subset_name = context.result['name']

    result = context.ptmscoutapp.post('http://localhost/experiments/28/subsets/share', {'saved-subset-select': context.subset_name})
    result.mustcontain(strings.share_subsets_page_title % (context.subset_name))

    m = re.search(r"http://localhost/experiments/28/subsets/share\?token=([a-zA-Z0-9]+)", result.body)
    context.share_link = m.group(0)

@then(u'other users should be able to view the shared subset')
def other_user_login_and_view(context):
    new_user = Bot(context.ptmscoutapp)
    new_user.register()
    new_user.activate()
    new_user.acquire_experiments([28])
    
    context.active_user.logout()
    new_user.login()
    
    result = context.ptmscoutapp.get(context.share_link)

    result.mustcontain(strings.share_subsets_page_title % (context.subset_name))
    result.mustcontain(strings.share_subset_success_message % ("Effects of HER2 overexpression on cell signaling networks governing proliferation and migration.", "http://localhost/experiments/28/subsets?annotation_set=%d" % (context.annotation_set_id)))

    result = new_user.fetch_subset_from_dataset_explorer(28, "%s: %s" % (context.active_user.username, context.subset_name), context.annotation_set_id)
    assertEqual( 26, result['foreground']['peptides'] )

@given(u'a user uploads a file containing non-mass spec experimental data')
def upload_dataset_for_annotation(context):
    context.active_user.login()
    
    filename = 'datasetExplorer_sample.txt'
    form = context.ptmscoutapp.get('/dataset/upload',status=200).form
    
    filename = os.path.join('tests','behave','data',filename)
    
    f = open(filename, 'rb')
    filecontents = f.read()
    
    form.set('load_type', "new")
    form.set('datasetname', "Test Data Set with 1 Error")
    form.set('datasetfile', (filename, filecontents))
    context.result = form.submit().follow()
    
    context.result = context.result.form.submit()
    if context.result.status != '302 Found':
        context.result.showbrowser()
    assert context.result.status == '302 Found'
    
    context.result = context.result.follow()
    
    context.form = context.result.form
    context.form.set('terms_of_use', "yes")
    
    
@then(u'the user should be sent an email with a link to the dataset which contains')
@patch('transaction.abort')
@patch('transaction.commit')
@patch('ptmscout.utils.mail.celery_send_mail')
def experiment_uploaded_check_email(context, patch_mail, patch_commit, patch_abort):
    from ptmscout.database import DBSession
    
    patch_commit.side_effect = bot.session_flush
    patch_abort.side_effect = bot.log_abort

    context.result = context.form.submit()

    exp_title="Test Data Set with 1 Error"
    DBSession.flush()

    my_experiments_page = "http://localhost/account/experiments"
    result = context.result

    result.mustcontain(strings.dataset_upload_started_page_title)
    result.mustcontain(strings.dataset_upload_started_message % (my_experiments_page))
    
    result = context.ptmscoutapp.get(my_experiments_page)
    
    result.mustcontain(exp_title)
    result.mustcontain("finished")

    m = re.search('http://localhost/experiments/([0-9]+)', str(result))
    if m == None:
        print 'http://localhost/experiments/[0-9]+ not found'
    exp_link = m.group(0)
    context.exp_id = int(m.group(1))

    exp = experiment.getExperimentById(context.exp_id, context.active_user.user)
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

@then(u'the user should be able to export their dataset with additional annotations')
@patch('transaction.abort')
@patch('transaction.commit')
@patch('ptmscout.utils.mail.celery_send_mail')
def download_annotated_experiment_and_check_fields(context, patch_mail, patch_commit, patch_abort):
    patch_commit.side_effect = bot.session_flush
    patch_abort.side_effect = bot.log_abort

    result = context.ptmscoutapp.get('http://localhost/experiments/%d/export?annotate=yes' % (context.exp_id))
    result.mustcontain(strings.experiment_export_started_page_title)


    context.mailargs = str(patch_mail.call_args_list)
    m = re.search(r'<a href="(.*)">here</a>', context.mailargs)

    download_url = m.group(1)
    result = context.ptmscoutapp.get(download_url)

    for row in context.table:
        field = row['field']
        num = int(row['elements'])
        
        cnt = 0
        cr = csv.DictReader(result.body.split("\n"), dialect='excel-tab')
        for crow in cr:
            crow[field] = crow[field].strip()

            if crow[field] != '':
                cnt += len(crow[field].split(";"))
        
        assert num == cnt, "expected %d but got %d, field %s is not as expected" % (num, cnt, field)


@then(u'the user should be able to view summary data')
def show_summaries(context):
    result = context.ptmscoutapp.get('http://localhost/experiments/%d/GO' % (context.exp_id))
    result = context.ptmscoutapp.get('http://localhost/experiments/%d/pfam' % (context.exp_id))
    result = context.ptmscoutapp.get('http://localhost/experiments/%d/predictions' % (context.exp_id))
    result = context.ptmscoutapp.get('http://localhost/experiments/%d/errors' % (context.exp_id))

@then(u'the user should be able to use the dataset explorer')
def show_explorer(context):
    result = context.ptmscoutapp.get('http://localhost/experiments/%d/subsets' % context.exp_id)
    
#    d = result.pyquery
#    json_base64 = d('#field-data').text()
#    response = json.loads( base64.b64decode(json_base64) )
    
    fq = [
          ['nop', ['quantitative', 'average:data:10', 'gt', '2']]
         ]
    
    result = context.active_user.query_dataset_explorer( context.exp_id, None, fq )
    assertEqual( 5, result['foreground']['peptides'] )
    
    context.result = result
    
@then(u'the user should not be able to use other experiment specific tools')
def user_forbidden_comparison_and_ambiguity(context):
    result = context.ptmscoutapp.get('http://localhost/experiments/%d' % (context.exp_id))
    assert result.body.find('Compare Datasets') == -1, "Compare Datasets should not be available"
    assert result.body.find('Report Ambiguity') == -1, "Report Ambiguity should not be available"
    
    context.ptmscoutapp.get('http://localhost/experiments/%d/compare' % (context.exp_id), status=403)
    context.ptmscoutapp.get('http://localhost/experiments/%d/ambiguity' % (context.exp_id), status=403)
   
@then(u'the user should be able to delete the dataset if needed')
def user_delete_dataset(context):
    result = context.ptmscoutapp.get('http://localhost/account/experiments')
    result = result.click(linkid='delete-dataset-0')
    result.mustcontain(strings.delete_experiment_confirm_message)

    result = result.forms['confirm'].submit()
    result.mustcontain(strings.delete_experiment_success_message)

    result = context.ptmscoutapp.get('http://localhost/experiments/%d' % (context.exp_id), status=403)


@given(u'a user uploads multiple clusterings for an experiment')
def user_upload_cluster_file(context):
    context.active_user.login()
    
    result = context.ptmscoutapp.get('/experiments/28/subsets')
    annote_file = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'experiment.28.test.clusters.tsv')
    result.forms[1].set('annotationname', "Annotation Set 1")
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

@when(u'the user chooses to export MCAM input files')
@patch('transaction.abort')
@patch('transaction.commit')
@patch('ptmscout.utils.mail.celery_send_mail')
def user_export_mcam(context, patch_mail, patch_commit, patch_abort):
    patch_commit.side_effect = bot.session_flush
    patch_abort.side_effect = bot.log_abort
    
    # submit the clusterings
    context.result.form.submit()

    from ptmscout.database import annotations
    from ptmscout.database import DBSession

    # submit the MCAM request
    result = context.ptmscoutapp.get('/experiments/28/mcam_enrichment')

    annot_set_id = int(result.form['annotationset'].options[1][0])
    print "User login:", context.active_user.user.id
    print "Selecting:", annot_set_id

    for annot_set in DBSession.query(annotations.AnnotationSet):
        print "Annotation set: (%d, '%s', %d)" % ( annot_set.id, annot_set.name, annot_set.experiment_id )
        print "Types:"
        for annot_type in annot_set.annotation_types:
            print "   ", annot_type.name, annot_type.type
        print "Permissions:"
        for perm in annot_set.permissions:
            print "Permission:", perm.user_id, perm.annotation_set_id

    annot_set = annotations.getUserAnnotations(annot_set_id, 28, context.active_user.user)
    print "Available:", "None" if annot_set is None else annot_set.name

    result.form['annotationset'] = result.form['annotationset'].options[1][0]
    result.form.set('correction', 'fdr')
    result.form.set('scansitecutoff', "3.0")
    result.form.set('domaincutoff', "1e-5")
    result.form.set('alpha', "0.05")
    result = result.form.submit()
    result.follow()

    context.mailargs = str(patch_mail.call_args_list)


@then(u'the user should be sent an email with a link to download their analysis')
def user_receive_email(context):
    assertContains(context.active_user.user.email, context.mailargs)
    assertContains(strings.annotation_upload_finished_subject, context.mailargs)
    assertContains("Annotations Loaded: 871", context.mailargs)
    assertContains("Errors Encountered: 1", context.mailargs)
    
    assertContains(strings.mcam_enrichment_finished_subject, context.mailargs)

@then(u'the user should be able to download an archive containing the MCAM files')
def user_download_mcam_files(context):
    m = re.search(r'<a href="(.*)">here</a>', context.mailargs)
    assert m != None, "Could not find link to MCAM results in: " + context.mailargs
    
    mcam_file = context.ptmscoutapp.get(m.group(1))
    
    headers = str(mcam_file.headers)
    m = re.search(r'filename="(.*)\.zip"', headers)
    assert m != None, "No filename found in header: " + headers
    
    zdir = m.group(1)
    zfilename = "%s.zip" % (zdir)
    
    tz = open(zdir + ".zip", 'wb')
    tz.write( mcam_file.body )
    tz.close()
    
    zf = zipfile.ZipFile(open(zfilename, 'r'))
    zf.extractall()
    
    numStructFilename = os.path.join(zdir, 'loadNumStruct.m')
    enrichBoolFilename = os.path.join(zdir, 'loadEnrichBool.m')
    numStructExpFilename = "tests/behave/data/datasetExplorer_numStruct_expected.txt"
    enrichBoolExpFilename = "tests/behave/data/datasetExplorer_enrichBool_expected.txt"

    if False:
        import shutil
        shutil.copyfile(numStructFilename, numStructExpFilename)
        shutil.copyfile(enrichBoolFilename, enrichBoolExpFilename)

    assertions.assertTextFilesEqual(numStructFilename, numStructExpFilename)
    assertions.assertTextFilesEqual(enrichBoolFilename, enrichBoolExpFilename)
    
    os.remove(numStructFilename)
    os.remove(enrichBoolFilename)
    os.removedirs(zdir)
    os.remove(zfilename)
