from behave import *
from ptmscout.database import upload
from ptmscout.config import strings

@given(u'a user is selecting experimental conditions for a dataset')
def open_experiment_conditions_page(context):
    context.active_user.login()
    
    form = context.ptmscoutapp.get("/upload", status=200).form
    result = context.active_user.load_datafile('datasetLoad_correctDataset.txt', form).follow()
    
    # submit the column settings form, as it should be filled out correctly
    result = result.form.submit(status=302).follow(status=200)
    
    # fill out the experimental metadata form to create an experiment
    result.form.set('experiment_name', 'some correct experiment')
    result.form.set('published', 'no')
    result.form.set('URL', '')
    result.form.set('description', 'a correct experiment for test purposes')
    result.form.set('ambiguous', 'no')
    result.form.set('notes','')
    
    context.result = result.form.submit().follow()
    
@then(u'give the user options to add cell, tissue, drug, stimulation, or environmental conditions')
def check_form_exists_with_cell_tissue_conditions(context):
    form = context.result.form
    form.set('0_type', 'tissue')
    form.set('0_value', 'BM-CD34')

    form.set('1_type', 'cell')
    form.set('1_value', 'HELA')
    
    form.set('2_type', 'drug')
    form.set('2_value', 'myostatin')
    
    form.set('3_type', 'stimulus')
    form.set('3_value', 'alot')
    
    form.set('4_type', 'environment')
    form.set('4_value', 'ph=5.5')
    
    result = form.submit()
    result.showbrowser()
    result = result.follow()
    result.mustcontain(strings.experiment_upload_confirm_message)

@then(u'automatically suggest suitable completions based on existing entries')
def check_autocomplete_function(context):
    field_names = set()
    for row in context.table:
        field_names.add(row['field_name'])
        
    auto_completions = {}
    
    for field in field_names:
        json = context.ptmscoutapp.get("/webservice/autocomplete/%s" % (field)).json
        auto_completions[field] = json[field]
        
    for row in context.table:
        field = row['field_name']
        suggestion = row['suggestion']
        
        assert suggestion in auto_completions[field]
    
    