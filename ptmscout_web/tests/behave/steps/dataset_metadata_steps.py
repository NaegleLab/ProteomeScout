from behave import *
from ptmscout.database import upload

def create_column(cnum, ctype, clabel=''):
    c = upload.SessionColumn()
    c.column_number = cnum
    c.type=ctype
    c.label=clabel
    return c

@given(u'a user is submitting publication information for a new dataset')
def setup_pre_existing_session(context):
    context.active_user.login()
    session = upload.Session()
    
    session.experiment_id = None
    session.load_type = 'new'
    session.parent_experiment = None
    session.change_description = ''
    session.data_file = 'some_dataset'
    session.user_id = context.active_user.user.id
    session.units = 'time(min)'
    session.stage = 'complete'
    session.save()
    
    context.result = context.ptmscoutapp.get('/upload/%d/metadata' % (session.id), status=200)

@when(u'the user submits a pubmed ID to fill out citation information')
def query_pubmed_for_citation(context):
    context.result.mustcontain("Published")
    context.result.mustcontain("Load Pubmed Citation")

    
#{'STAT': 'MEDLINE', 'IP': '3', 'JT': 'Briefings in bioinformatics', 'DA': '20020916', 'FAU': ['Mangalam, Harry'], 'DP': '2002 Sep', 'OWN': 'NLM', 'PT': ['Journal Article'], 'LA': ['eng'], 'CRDT': ['2002/09/17 10:00'], 'DCOM': '20030606', 'LR': '20041117', 'PG': '296-302', 'TI': 'The Bio* toolkits--a brief overview.', 'PL': 'England', 'TA': 'Brief Bioinform', 'JID': '100912837', 'AB': 'Bioinformatics research is often difficult to do with commercial software. The Open Source BioPerl, BioPython and Biojava projects provide toolkits with multiple functionality that make it easier to create customised pipelines or analysis. This review briefly compares the quirks of the underlying languages and the functionality, documentation, utility and relative advantages of the Bio counterparts, particularly from the point of view of the beginning biologist programmer.', 'AD': 'tacg Informatics, Irvine, CA 92612, USA. hjm@tacgi.com', 'VI': '3', 'IS': '1467-5463 (Print) 1467-5463 (Linking)', 'AU': ['Mangalam H'], 'MHDA': '2003/06/07 05:00', 'MH': ['*Computational Biology', 'Computer Systems', 'Humans', 'Internet', '*Programming Languages', '*Software', 'User-Computer Interface'], 'EDAT': '2002/09/17 10:00', 'SO': 'Brief Bioinform. 2002 Sep;3(3):296-302.', 'SB': 'IM', 'PMID': '12230038', 'PST': 'ppublish'}
@then(u'automatically fill out the publication information from the pubmed record')
def check_pubmed_result(context):
    context.result = context.ptmscoutapp.get('/webservice/pubmed/12230038', status=200)
       
    cite_dict = context.result.json
    
    cite_dict['authors'] = "Mangalam H"
    cite_dict['journal'] = "Briefings in bioinformatics"
    cite_dict['publication_month'] = "september"
    cite_dict['publication_year'] = "2002"
    cite_dict['volume'] = "3"
    cite_dict['page_start'] = "296"
    cite_dict['page_end'] = "302"
    
