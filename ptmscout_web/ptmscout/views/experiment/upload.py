from pyramid.view import view_config
from ptmscout.views.layout import site_layout
from ptmscout.config import strings
import base64
import json
from ptmscout.utils import webutils
from pyramid.httpexceptions import HTTPFound

def create_experiment_and_start_upload(field_dict, tmp_file):
    pass

def verify_datafile(request):
    pass


def get_enum_fields(request, users_experiments):
    return {'parent_experiment': set([ str(e.id) for e in users_experiments ]),
            'published': set(["no","yes"]),
            'ambiguous': set(["no","yes"]),
            'publication_month':set(["january","february", "march", "april",
                                 "may","june","july","august","september",
                                 "october","november","december"])}

def get_numeric_fields(request):
    return ['pmid',
            'publication_year',
            'volume']

def get_required_fields(request):
    req_fields = ['load_type', 
                  'data_file', 
                  'description', 
                  'experiment_name', 
                  'published', 
                  'author_contact'
                  'ambiguous']
    
    if (webutils.post(request, 'load_type', "") != "new"):
        req_fields.append('parent_experiment')
    
    if (webutils.post(request, 'load_type', "") == "extension"):
        req_fields.append('change_description')
        
    if (webutils.post(request, 'published', "") == "yes"):
        req_fields.append('authors')
        req_fields.append('journal')
        req_fields.append('publication_year')
        req_fields.append('publication_month')
        req_fields.append('volume')
        req_fields.append('page_start')
        req_fields.append('page_end')
    return req_fields

def check_required_fields(request, users_experiments):
    req_fields = get_required_fields(request)
    
    field_name_dict={'pmid':"PubMed ID", 
                     'publication_year': "Publication Year",
                     'publication_month': "Publication Month", 
                     'parent_experiment': "Parent Experiment", 
                     'published': "Published",
                     'ambiguous': "Possibly ambiguous accessions",
                     'author_contact' : "Author Contact Email"}
    
    field_dict = {}
    for field in request.POST:
        field_dict[field] = webutils.get(request, field, "").strip()

    if('terms_of_use' not in field_dict):
        return False, strings.failure_reason_terms_of_use_not_accepted, field_dict
    
    del field_dict['terms_of_use']
    
    for field in req_fields:
        if field_dict[field] == "":    
            return False, strings.failure_reason_required_fields_cannot_be_empty, field_dict
        
    numeric_fields = get_numeric_fields(request)
    
    for field in numeric_fields:
        if field_dict[field] not in req_fields and field_dict[field] == "":
            continue
        try:
            int(field_dict[field])
        except:
            return False, strings.failure_reason_field_must_be_numeric % field_name_dict[field], field_dict
    
    enum_fields = get_enum_fields(request, users_experiments)
    
    for field in enum_fields:
        if field_dict[field] not in req_fields and field_dict[field] == "":
            continue
        if not field_dict[field] in enum_fields[field]:
            return False, strings.failure_reason_field_value_not_valid % field_name_dict[field], field_dict
    
    return True, None, field_dict
    

def experiment_to_dict(exp):
    expd = {}
    
    for key in exp.__dict__:
        if key[0] != "_":
            expd[key] = exp.__dict__[key]
    
    return expd

@view_config(route_name='upload', renderer='ptmscout:templates/accounts/upload.pt', permission='private')
def user_upload(request):
    submitted = webutils.post(request,'submitted',"false")
        
    users_experiments = [ p.experiment for p in request.user.permissions if p.access_level=='owner' ]
    dict_exps = [ experiment_to_dict(exp) for exp in users_experiments ]
    reason = None
    field_dict = {}
    
    if(submitted == "true"):
        passed, error, field_dict = check_required_fields(request, users_experiments)
        
        if not passed:
            reason = error
        else:
            passed, result = verify_datafile(request)
            
            if not passed:
                reason = result
            else:
                exp_id = create_experiment_and_start_upload(field_dict, result)
                return HTTPFound(request.application_url + "/upload/" + str(exp_id))
    
    return {'pageTitle': strings.upload_page_title,
            'user_experiments': dict_exps,
            'json_user_data': base64.b64encode(json.dumps(dict_exps)),
            'reason':reason,
            'formfields':base64.b64encode(json.dumps(field_dict))}
