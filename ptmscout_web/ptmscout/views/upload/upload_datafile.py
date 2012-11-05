from pyramid.view import view_config
import time
from ptmscout.config import settings, strings
import os
from ptmscout.utils import webutils, forms
from pyramid.httpexceptions import HTTPFound, HTTPForbidden
from ptmscout.database import upload
import logging
log = logging.getLogger(__name__)

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
    
    os.system("mac2unix -q %s" % os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, exp_file))
    os.system("dos2unix -q %s" % os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, exp_file))

    return exp_file

def create_schema(request, users_experiments):
    schema = forms.FormSchema()
    
    parent_experiment_options = [(str(e.id), e.name) for e in users_experiments]
    
    schema.add_radio_field('load_type', "Load Type", [('new',"New"),('append',"Append"),('reload',"Reload"),('extension',"Extension")])
    schema.add_select_field('parent_experiment', 'Parent Experiment', parent_experiment_options)
    schema.add_textarea_field('change_description', "Change Description", 43, 5)
    schema.add_file_upload_field('data_file', 'Input Data File')
    
    schema.set_field_required_condition('change_description', 'load_type', lambda pval: pval == "extension")
    schema.set_field_required_condition('parent_experiment', 'load_type', lambda pval: pval != "new")
    schema.set_required_field('load_type')
    schema.set_required_field('data_file')
    
    schema.parse_fields(request)
    
    return schema

    
@view_config(route_name='upload', renderer='ptmscout:/templates/upload/upload_datafile.pt', permission='private')
def upload_data_file(request):
    submitted = webutils.post(request, 'submitted', "false") == "true"
    users_experiments = [ p.experiment for p in request.user.permissions if p.access_level=='owner' ]
        
    errors = []
    schema = create_schema(request, users_experiments)
    
    if submitted:
        errors = forms.FormValidator(schema).validate()
        
        if len(errors) == 0:
            output_file = save_data_file(request)
            session_id = create_session(request, output_file)
            return HTTPFound(request.application_url + "/upload/%d/config" % (session_id))

    return {'pageTitle': strings.upload_page_title,
            'header': strings.upload_page_header,
            'formrenderer': forms.FormRenderer(schema),
            'errors':errors}