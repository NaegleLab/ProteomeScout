from behave import *
import json
from ptmscout.database import experiment, upload
from ptmscout.config import strings
from tests.behave.steps.assertions import assertEqual

# ctype in set(['data','stddev','accession','peptide','species','modification','run', 'none'])
def create_column(cnum, ctype, clabel=''):
    c = upload.SessionColumn()
    c.column_number = cnum
    c.type=ctype
    c.label=clabel
    return c

def setup_pre_existing_experiment_and_session(context):
    session = upload.Session()
    exp = experiment.Experiment()
    
    exp.name = "Some experiment"
    exp.author = "Some author"
    exp.dataset = "some_dataset"
    exp.description = "Some description"
    exp.errorLog = "Some log"
    exp.status = "loaded"
    exp.submitter_id = context.active_user.user.id
    exp.published = 0
    exp.public = 0
    exp.ambiguity = 0
    exp.export = 1
    exp.grantPermission(context.active_user.user, 'owner')
    exp.saveExperiment()
    
    session.experiment_id = exp.id
    session.load_type = 'new'
    session.parent_experiment = None
    session.change_description = ''
    session.data_file = 'some_dataset'
    session.user_id = context.active_user.user.id
    session.units = 'time(min)'
    session.stage = 'complete'
    
    c1 = create_column(0, 'accession')
    c2 = create_column(2, 'modification')
    c3 = create_column(3, 'peptide')
    c4 = create_column(4, 'run')
    
    c5 = create_column(5, 'data', '0')
    c6 = create_column(6, 'data', '30')
    c7 = create_column(13, 'stddev', '0')
    c8 = create_column(14, 'stddev', '30')
    
    session.columns.extend([c1,c2,c3,c4,c5,c6,c7,c8])
    session.save()
    
    return exp.id


@given(u'a user is loading a dataset')
def access_data_upload_form(context):
    context.active_user.login()
    context.preexisting_exp = setup_pre_existing_experiment_and_session(context)
    context.result = context.ptmscoutapp.get('/upload', status=200)
    context.form = context.result.form

    
@when(u'some amino acids in lower case do not match any species possibilities for that type of modification')
def load_data_file_with_bad_mod_for_amino_acid(context):
    context.result = context.active_user.load_datafile('datasetLoad_badAminoAcid.txt', context.form).follow()

@when(u'identical peptide/protein pairs exist with different data and no explicit run column')
def load_data_file_with_no_run_column(context):
    context.result = context.active_user.load_datafile('datasetUI_replicateEntriesWithNoRunCol.txt', context.form).follow()

@when(u'the dataset exists and new entries are being appended')
def load_data_file_when_appending(context):
    context.result = context.active_user.load_datafile('datasetLoad_correctDatasetNoHeader.txt', context.form, load_type='append', parent_experiment=str(context.preexisting_exp)).follow()

@when(u'the user wants the dataset to replace another dataset')
def load_data_file_when_replacing(context):
    context.result = context.active_user.load_datafile('datasetLoad_correctDatasetNoHeader.txt', context.form, load_type='reload', parent_experiment=str(context.preexisting_exp)).follow()

@when(u'the user wants the dataset to extend another dataset')
def load_data_file_when_extending(context):
    context.result = context.active_user.load_datafile('datasetLoad_correctDatasetNoHeader.txt', context.form, load_type='extension', parent_experiment=str(context.preexisting_exp), change_description="Some set of changes are happening").follow()


@when(u'non-alphabetic characters appear in some of the entries')
def load_data_with_non_alphabetics_in_peptide(context):
    context.result = context.active_user.load_datafile('datasetUI_badPep.txt', context.form).follow()

@when(u'data columns exist with declarations acc, MOD_TYPE, pep, data:type:value and matching stddev:type:value')
def load_data_when_good_headers(context):
    context.result = context.active_user.load_datafile('datasetUI_explicitDataHeaders.txt', context.form).follow()


@when(u'headers cannot be identified that conform to acc, MOD_TYPE, data, stddev and pep')
def load_data_when_unrecognized_headers(context):
    context.result = context.active_user.load_datafile('datasetUI_unrecognizedHeaders.txt', context.form).follow()


@when(u'the user submits a pubmed ID to fill out citation information')
def query_pubmed_for_citation(context):
    context.result = context.ptmscoutapp.get('/webservice/pubmed/12230038', status=200)


def assert_column_type(form, i, t, label=None):
    form_elem = form.get('column_%d_type' % (i))
    print form_elem.options
    selected_options = [ name for (name, selected) in form_elem.options if selected ]
    
    assertEqual([t], selected_options)
    if label != None: 
        assertEqual(label, form.get('column_%d_label' % (i)).value)

@then(u'show the user column assignments and data column assignments with type and value fields')
def check_fields_prepopulated_properly(context):
    assert_column_type(context.result.form, 0, 'none')
    assert_column_type(context.result.form, 1, 'accession')
    assert_column_type(context.result.form, 2, 'none')
    assert_column_type(context.result.form, 3, 'modification')
    assert_column_type(context.result.form, 4, 'none')
    assert_column_type(context.result.form, 5, 'none')
    assert_column_type(context.result.form, 6, 'none')
    assert_column_type(context.result.form, 7, 'none')
    assert_column_type(context.result.form, 8, 'none')
    assert_column_type(context.result.form, 9, 'peptide')
    
    assert_column_type(context.result.form, 10, 'data', '0')
    assert_column_type(context.result.form, 11, 'data', '1')
    assert_column_type(context.result.form, 12, 'data', '5')
    
    assert_column_type(context.result.form, 13, 'stddev', '1')
    assert_column_type(context.result.form, 14, 'stddev', '5')
    
    assertEqual('time(min)', context.result.form.get('units').value)

@then(u'show the user their headers')
def check_fields_not_prepopulated_correct_header_titles_shown(context):
    headers = "NAME    record id    Gene Symbol    modification    RSD    SITE_GRP_ID    SPECIES    MW (kD)    IN_DOMAIN    tryps    measurement 1    measurement 2    measurement 3".split("    ")
    
    for header in headers:
        context.result.mustcontain(header)
        

@then(u'show the user that incorrect modifications types have been detected')
def check_errors_reported_bad_modification_for_amino_acid(context):
    result = context.result.form.submit()
    result.mustcontain(strings.experiment_upload_warning_modifications_do_not_match_amino_acids)

@then(u'show the user that replicate data appears to have been detected')
def check_errors_reported_replicate_data(context):
    result = context.result.form.submit()
    result.mustcontain(strings.experiment_upload_warning_no_run_column)
    
@then(u'show the user that bad peptide strings have been detected')
def check_errors_reported_bad_peptide_definition(context):
    result = context.result.form.submit()
    result.mustcontain(strings.experiment_upload_warning_peptide_column_contains_bad_peptide_strings)

@then(u'pre-populate all fields of data loading with assignments from original dataset')
def check_fields_prepopulated_properly_from_session_data(context):
    assert_column_type(context.result.form, 0, 'accession') 
    assert_column_type(context.result.form, 2, 'modification') 
    assert_column_type(context.result.form, 3, 'peptide') 
    assert_column_type(context.result.form, 4, 'run')
    
    assert_column_type(context.result.form, 5, 'data') 
    assert_column_type(context.result.form, 6, 'data')
    assert_column_type(context.result.form, 5, 'data', '0')
    assert_column_type(context.result.form, 6, 'data', '30')
    
    assert_column_type(context.result.form, 13, 'stddev')
    assert_column_type(context.result.form, 14, 'stddev')
    assert_column_type(context.result.form, 13, 'stddev', '0')
    assert_column_type(context.result.form, 14, 'stddev', '30')
    
    assertEqual('time(min)', context.result.form.get('units').value)

    
#{'STAT': 'MEDLINE', 'IP': '3', 'JT': 'Briefings in bioinformatics', 'DA': '20020916', 'FAU': ['Mangalam, Harry'], 'DP': '2002 Sep', 'OWN': 'NLM', 'PT': ['Journal Article'], 'LA': ['eng'], 'CRDT': ['2002/09/17 10:00'], 'DCOM': '20030606', 'LR': '20041117', 'PG': '296-302', 'TI': 'The Bio* toolkits--a brief overview.', 'PL': 'England', 'TA': 'Brief Bioinform', 'JID': '100912837', 'AB': 'Bioinformatics research is often difficult to do with commercial software. The Open Source BioPerl, BioPython and Biojava projects provide toolkits with multiple functionality that make it easier to create customised pipelines or analysis. This review briefly compares the quirks of the underlying languages and the functionality, documentation, utility and relative advantages of the Bio counterparts, particularly from the point of view of the beginning biologist programmer.', 'AD': 'tacg Informatics, Irvine, CA 92612, USA. hjm@tacgi.com', 'VI': '3', 'IS': '1467-5463 (Print) 1467-5463 (Linking)', 'AU': ['Mangalam H'], 'MHDA': '2003/06/07 05:00', 'MH': ['*Computational Biology', 'Computer Systems', 'Humans', 'Internet', '*Programming Languages', '*Software', 'User-Computer Interface'], 'EDAT': '2002/09/17 10:00', 'SO': 'Brief Bioinform. 2002 Sep;3(3):296-302.', 'SB': 'IM', 'PMID': '12230038', 'PST': 'ppublish'}
@then(u'automatically fill out the publication information from the pubmed record')
def check_pubmed_result(context):
    cite_dict = json.loads(context.result)
    
    cite_dict['authors'] = "Mangalam H"
    cite_dict['journal'] = "Briefings in bioinformatics"
    cite_dict['publication_month'] = "september"
    cite_dict['publication_year'] = "2002"
    cite_dict['volume'] = "3"
    cite_dict['page_start'] = "296"
    cite_dict['page_end'] = "302"
    
