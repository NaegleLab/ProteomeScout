from pyramid.view import view_config
from ptmscout.config import strings
import base64
import json
from ptmscout.utils import webutils
from pyramid.httpexceptions import HTTPFound
from ptmscout.database import experiment, upload

def create_experiment(field_dict, session, current_user):
    if(session.load_type == 'reload' or session.load_type == 'append'):
        nexp = experiment.getExperimentById(session.parent_experiment, current_user)
    else: 
        nexp = experiment.Experiment()
        nexp.public = 0
        nexp.experiment_id = None

    nexp.name = field_dict['experiment_name']
    nexp.description = field_dict['description']
    
    nexp.dataset = session.data_file
    nexp.status='preload'
    
    nexp.published = 1 if field_dict['published'] == "yes" else 0
    nexp.ambiguity = 1 if field_dict['ambiguous'] == "yes" else 0
    nexp.URL = field_dict['URL']
    nexp.export = 1

    nexp.contact = None
    nexp.author = None
    nexp.journal = None
    nexp.volume = None
    nexp.page_end = None
    nexp.page_start = None
    nexp.publication_month = None
    nexp.publication_year = None
    nexp.errorLog = None
    nexp.PMID = None
    
    nexp.submitter_id = current_user.id
    
    if(session.load_type == 'extension'):
        nexp.experiment_id = session.parent_experiment
        
    if(field_dict['published'] == "yes"):
        nexp.contact = field_dict['author_contact']
        nexp.author = field_dict['authors']
        
        nexp.volume = int(field_dict['volume'])
        nexp.page_start = field_dict['page_start']
        nexp.page_end = field_dict['page_end']
        
        nexp.publication_month = field_dict['publication_month']
        nexp.publication_year = int(field_dict['publication_year'])
        
        nexp.journal = field_dict['journal']
        if(field_dict['pmid'] != ""):
            nexp.PMID = int(field_dict['pmid'])
    
    
    nexp.grantPermission(current_user, 'owner')
    nexp.saveExperiment()
    
    return nexp
    
def create_experiment_and_mark_status(field_dict, session_id, current_user):
    session = upload.getSessionById(session_id, current_user)
    nexp = create_experiment(field_dict, session, current_user)
    
    session.experiment_id = nexp.id
    session.stage = 'confirm'
    session.save()

def get_enum_fields(request):
    return {'published': set(["no","yes"]),
            'ambiguous': set(["no","yes"]),
            'publication_month':set(["january","february", "march", "april",
                                 "may","june","july","august","september",
                                 "october","november","december"])}

def get_numeric_fields(request):
    return ['pmid',
            'publication_year',
            'volume']

def get_required_fields(request):
    req_fields = ['description', 
                  'experiment_name', 
                  'published', 
                  'ambiguous']
    
    if (webutils.post(request, 'published', "") == "yes"):
        req_fields.append('author_contact')
        req_fields.append('authors')
        req_fields.append('journal')
        req_fields.append('publication_year')
        req_fields.append('publication_month')
        req_fields.append('volume')
        req_fields.append('page_start')
        req_fields.append('page_end')
    return req_fields

def check_required_fields(request):
    req_fields = get_required_fields(request)
    
    field_name_dict={'pmid':"PubMed ID", 
                     'publication_year': "Publication Year",
                     'publication_month': "Publication Month", 
                     'published': "Published",
                     'ambiguous': "Possibly ambiguous accessions",
                     'author_contact' : "Author Contact Email",
                     'page_start': "Page Start",
                     'page_end': "Page End",
                     'author_contact':"Author Contact",
                     'authors': "Authors",
                     'journal': "Journal",
                     'volume': "Volume",
                     'description': "Description",
                     'experiment_name': "Experiment Name"}
    
    field_dict = {}
    for field in request.POST:
        field_dict[field] = webutils.post(request, field, "")
        if(isinstance(field_dict[field], str)):
            field_dict[field] = field_dict[field].strip()

    if('terms_of_use' not in field_dict):
        return False, strings.failure_reason_terms_of_use_not_accepted, field_dict
    
    del field_dict['terms_of_use']
    
    for field in req_fields:
        if field not in field_dict or field_dict[field] == "":
            return False, strings.failure_reason_required_fields_cannot_be_empty % field_name_dict[field], field_dict
        
    numeric_fields = get_numeric_fields(request)
    
    for field in numeric_fields:
        if field_dict[field] not in req_fields and field_dict[field] == "":
            continue
        try:
            int(field_dict[field])
        except:
            return False, strings.failure_reason_field_must_be_numeric % field_name_dict[field], field_dict
    
    enum_fields = get_enum_fields(request)
    
    for field in enum_fields:
        if field not in req_fields and field_dict[field] == "":
            continue
        if not field_dict[field] in enum_fields[field]:
            return False, strings.failure_reason_field_value_not_valid % field_name_dict[field], field_dict
    
    return True, None, field_dict
    
@view_config(route_name='upload_metadata', renderer='ptmscout:templates/upload/upload_metadata.pt', permission='private')
def upload_metadata(request):
    submitted = webutils.post(request,'submitted',"false")
    session_id = int(request.matchdict['id'])
    
    reason = None
    field_dict = {}
    
    if(submitted == "true"):
        passed, error, field_dict = check_required_fields(request)
        
        if not passed:
            reason = error
        else:
            # need to get exp_file from upload session records
            create_experiment_and_mark_status(field_dict, session_id, request.user)
            return HTTPFound(request.application_url + "/upload/%d/confirm" % session_id)
    
    if 'data_file' in field_dict:
        del field_dict['data_file']
    
    return {'pageTitle': strings.upload_page_title,
            'reason':reason,
            'formfields':base64.b64encode(json.dumps(field_dict))}
