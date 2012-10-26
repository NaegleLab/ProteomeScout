from pyramid.view import view_config
import time
from ptmscout.config import settings, strings
import os
from ptmscout.utils import webutils
import base64
import json
from pyramid.httpexceptions import HTTPFound, HTTPForbidden
from ptmscout.database import upload

def create_session(request, exp_file):
    session = upload.Session()
    
    session.user_id = request.user.id
    session.data_file = exp_file
    session.load_type = request.POST['load_type'].strip()
    session.parent_experiment = None
    session.stage = 'config'
    session.change_description = ''
    
    if session.load_type != 'new':
        session.parent_experiment = int(request.POST['parent_experiment'])
    
    if session.load_type == 'extension':
        session.change_description = request.POST['change_description']
    
    session.save()
    
    return session.id
    
def verify_data_file(exp_file):
    ifile = open(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, exp_file), 'rb')
    
    header = ifile.readline()
    
    return None if header.count("\t") >= 3 else strings.failure_reason_experiment_file_not_enough_columns 
    

def save_data_file(request):
    exp_file = "experiment_data" + str(time.time())
    
    input_file = request.POST['data_file'].file
    
    output_file = open(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, exp_file), 'wb')
    
    input_file.seek(0)
    while 1:
        data = input_file.read(2<<16)
        if not data:
            break
        output_file.write(data)
    output_file.close()
    
    os.system("mac2unix -q %s" % os.path.join(settings.experiment_data_file_path, exp_file))
    os.system("dos2unix -q %s" % os.path.join(settings.experiment_data_file_path, exp_file))
    
    error = verify_data_file(exp_file)
    
    if error == None:
        return True, exp_file
    
    return False, error


def experiment_to_dict(exp):
    expd = {}
    
    for key in exp.__dict__:
        if key[0] != "_":
            expd[key] = exp.__dict__[key]
    
    return expd

def get_required_fields(request):
    req_fields = ['data_file']
    if (webutils.post(request, 'load_type', "") != "new"):
        req_fields.append('parent_experiment')
    
    if (webutils.post(request, 'load_type', "") == "extension"):
        req_fields.append('change_description')
        
    return req_fields

def check_required_fields(request, users_experiments):
    field_name_dict = {'parent_experiment': "Parent Experiment",
                       'load_type': "Load Type",
                       'change_description': "Change Description"}
    
    field_dict = {}
    for field in request.POST:
        field_dict[field] = webutils.post(request, field, "")
        if(isinstance(field_dict[field], str)):
            field_dict[field] = field_dict[field].strip()
    
    req_fields = get_required_fields(request)
    
    for field in req_fields:
        if field not in field_dict or field_dict[field] == "":
            return False, strings.failure_reason_required_fields_cannot_be_empty % field_name_dict[field], field_dict
        
    valid_load_types = set(['new','append','reload','extension'])
    experiment_ids = set([str(e.id) for e in users_experiments])
    if field_dict['load_type'] not in valid_load_types or \
        (field_dict['load_type'] != 'new' and field_dict['parent_experiment'] not in experiment_ids):
        return False, strings.failure_reason_field_value_not_valid % field_name_dict['parent_experiment'], field_dict
    
    return True, None, field_dict
        
    
@view_config(route_name='upload', renderer='ptmscout:/templates/upload/upload.pt')
def upload_data_file(request):
    submitted = webutils.post(request, 'submitted', "false")
    
    if(request.user == None):
        raise HTTPForbidden()
    
    users_experiments = [ p.experiment for p in request.user.permissions if p.access_level=='owner' ]    
    dict_exps = [ experiment_to_dict(exp) for exp in users_experiments ]
    reason = None
    form_fields = {}
    
    if submitted == "true":
        success, reason, form_fields = check_required_fields(request, users_experiments)
        
        if success:
            success, result = save_data_file(request)
            
            if success:
                session_id = create_session(request, result)
                return HTTPFound(request.application_url + "/upload/%d/config" % (session_id))
            else:
                reason = result
        
    
    return {'pageTitle': strings.upload_page_title,
            'header': strings.upload_page_header,
            'user_experiments': dict_exps,
            'formfields': base64.b64encode(json.dumps(form_fields)),
            'reason':reason}