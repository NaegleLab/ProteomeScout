from pyramid.view import view_config
from ptmscout.config import strings, settings
import base64
import json
from ptmscout.utils import webutils
from pyramid.httpexceptions import HTTPFound
import time
import os
import datetime
from ptmscout.database import experiment

def create_experiment_and_start_upload(field_dict, exp_file, current_user):
    if(field_dict['load_type'] == 'reload'):
        nexp = experiment.getExperimentById(int(field_dict['parent_experiment']), current_user)
    else: 
        nexp = experiment.Experiment()
        nexp.public = 0
        nexp.experiment_id = None    

    nexp.name = field_dict['experiment_name']
    nexp.description = field_dict['description']
    
    nexp.dataset = exp_file
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
    
    if(field_dict['load_type'] == 'extension'):
        nexp.experiment_id = int(field_dict['parent_experiment'])
        
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
    return nexp.id

def verify_datafile(exp_file):
    infile = open(os.path.join(settings.experiment_data_file_path, exp_file), 'rb')
    header = [ col.strip().lower() for col in infile.readline().strip().split("\t") ]
    infile.close()
    
    acc_columns = []
    pep_columns = []
    
    for col in header:
        if col == "acc" or col.find("acc:") == 0:
            acc_columns.append(col)
        if col == "peptide" or col.find("pep:") == 0:
            pep_columns.append(col)
    
    if len(acc_columns) == 0:
        return strings.failure_reason_experiment_header_no_acc_column
    if len(acc_columns) > 1:
        return strings.failure_reason_experiment_header_multiple_acc_columns
    
    if len(pep_columns) == 0:
        return strings.failure_reason_experiment_header_no_peptide_column
    if len(pep_columns) > 1:
        return strings.failure_reason_experiment_header_multiple_peptide_column
    
    return None

def save_datafile(request):
    exp_file = "experiment_data" + str(time.time())
    
    input_file = request.POST['data_file'].file
    output_file = open(os.path.join(settings.experiment_data_file_path, exp_file), 'wb')
    
    input_file.seek(0)
    while 1:
        data = input_file.read(2<<16)
        if not data:
            break
        output_file.write(data)
    output_file.close()
    
    os.system("mac2unix -q %s" % os.path.join(settings.experiment_data_file_path, exp_file))
    os.system("dos2unix -q %s" % os.path.join(settings.experiment_data_file_path, exp_file))
    
    error = verify_datafile(exp_file)
    
    if error == None:
        return True, exp_file
    else:
        return False, error


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
                  'ambiguous']
    
    if (webutils.post(request, 'load_type', "") != "new"):
        req_fields.append('parent_experiment')
    
    if (webutils.post(request, 'load_type', "") == "extension"):
        req_fields.append('change_description')
        
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

def check_required_fields(request, users_experiments):
    req_fields = get_required_fields(request)
    
    field_name_dict={'pmid':"PubMed ID", 
                     'publication_year': "Publication Year",
                     'publication_month': "Publication Month", 
                     'parent_experiment': "Parent Experiment", 
                     'published': "Published",
                     'ambiguous': "Possibly ambiguous accessions",
                     'author_contact' : "Author Contact Email",
                     'page_start': "Page Start",
                     'page_end': "Page End",
                     'author_contact':"Author Contact",
                     'authors': "Authors",
                     'journal': "Journal",
                     'volume': "Volume",
                     'change_description': "Change Description",
                     'load_type': "Load Type",
                     'data_file': "Data File",
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
        if field_dict[field] == "":
            return False, strings.failure_reason_required_fields_cannot_be_empty % field_name_dict[field], field_dict
        
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
            passed, result = save_datafile(request)
            
            if not passed:
                reason = result
            else:
                exp_id = create_experiment_and_start_upload(field_dict, result, request.user)
                return HTTPFound(request.application_url + "/upload/" + str(exp_id))
    
    if 'data_file' in field_dict:
        del field_dict['data_file']
    
    dthandler = lambda obj: obj.isoformat() if isinstance(obj, datetime.datetime) else obj
    return {'pageTitle': strings.upload_page_title,
            'user_experiments': dict_exps,
            'json_user_data': base64.b64encode(json.dumps(dict_exps, default=dthandler)),
            'reason':reason,
            'formfields':base64.b64encode(json.dumps(field_dict))}
