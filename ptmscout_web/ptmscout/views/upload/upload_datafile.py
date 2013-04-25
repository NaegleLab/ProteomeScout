from pyramid.view import view_config
from ptmscout.config import strings, settings
from ptmscout.utils import webutils, forms, uploadutils
from pyramid.httpexceptions import HTTPFound
from ptmscout.database import upload
import logging
log = logging.getLogger(__name__)

def create_session(request, exp_file):
    session = upload.Session()
    
    session.user_id = request.user.id
    session.data_file = exp_file
    session.resource_type = 'experiment'
    session.load_type = request.POST['load_type'].strip()
    session.parent_experiment = None
    session.stage = 'config'
    session.change_name = ''
    session.change_description = ''
    
    if session.load_type != 'new':
        session.parent_experiment = int(request.POST['parent_experiment'])
    
    if session.load_type == 'extension':
        session.change_name = request.POST['change_name']
        session.change_description = request.POST['change_description']
    
    session.save()
    
    return session.id
    
def create_schema(request, users_experiments):
    schema = forms.FormSchema()
    
    parent_experiment_options = [(str(e.id), e.name) for e in users_experiments]
    
    schema.add_radio_field('load_type', "Load Type", [('new',"New"),('append',"Append"),('reload',"Reload"),('extension',"Extension")])

    schema.add_select_field('parent_experiment', 'Parent Experiment', parent_experiment_options)
    schema.add_text_field('change_name', "Extension Title", width=55)
    schema.add_textarea_field('change_description', "Description of Extension", 43, 5)
    schema.add_file_upload_field('data_file', 'Input Data File')
    
    schema.set_field_required_condition('change_name', 'load_type', lambda pval: pval == "extension")
    schema.set_field_required_condition('change_description', 'load_type', lambda pval: pval == "extension")
    schema.set_field_required_condition('parent_experiment', 'load_type', lambda pval: pval != "new")
    schema.set_required_field('load_type')
    schema.set_required_field('data_file')
    
    schema.parse_fields(request)
    
    return schema

    
@view_config(route_name='upload', renderer='ptmscout:/templates/upload/upload_datafile.pt', permission='private')
def upload_data_file(request):
    submitted = webutils.post(request, 'submitted', "false") == "true"
    users_experiments = [ e for e in request.user.myExperiments() if e.ready() ]
    
    errors = []
    schema = create_schema(request, users_experiments)
    
    if submitted:
        errors = forms.FormValidator(schema).validate()
        
        if len(errors) == 0:
            output_file = uploadutils.save_data_file(request.POST['data_file'], settings.experiment_files_prefix)
            session_id = create_session(request, output_file)
            return HTTPFound(request.application_url + "/upload/%d/config" % (session_id))

    return {'pageTitle': strings.upload_page_title,
            'header': strings.upload_page_header,
            'formrenderer': forms.FormRenderer(schema),
            'errors':errors}
