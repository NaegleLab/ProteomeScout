from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound
from ptmscout.config import settings, strings
from ptmscout.utils import uploadutils, forms
from ptmscout.database import experiment, upload

def create_session(experiment_id, annotation_filename, user):
    session = upload.Session()
    session.data_file = annotation_filename
    session.experiment_id = experiment_id
    
    session.load_type = 'annotations'
    session.stage = 'config'
    session.user_id = user.id
    session.change_description = ""
    
    session.save()

    return session.id

def create_schema(request):
    schema = forms.FormSchema()
    
    schema.add_file_upload_field('annotationfile', 'Input Data File')
    schema.set_required_field('annotationfile')
    schema.parse_fields(request)
    
    return schema

@view_config(route_name='experiment_annotate', request_method='POST',  renderer='ptmscout:/templates/experiments/experiment_annotate.pt', permission='private')
def upload_annotation_file_POST(request):
    experiment_id = int(request.matchdict['id'])
    exp = experiment.getExperimentById(experiment_id, request.user)
    
    schema = create_schema(request)
    errors = forms.FormValidator(schema).validate()
    
    if len(errors) == 0:
        output_file = uploadutils.save_data_file(request.POST['annotationfile'], settings.annotation_files_prefix)
        job_id = create_session(experiment_id, output_file, request.user)
        return HTTPFound(request.route_url('configure_annotations', id=experiment_id, sid=job_id))

    return {'experiment': exp,
            'pageTitle': strings.upload_annotations_page_title,
            'formrenderer': forms.FormRenderer(schema),
            'errors':errors}
    
@view_config(route_name='experiment_annotate', request_method='GET',  renderer='ptmscout:/templates/experiments/experiment_annotate.pt')
def upload_annotation_file_GET(request):
    experiment_id = int(request.matchdict['id'])
    exp = experiment.getExperimentById(experiment_id, request.user)
    
    schema = create_schema(request)
    
    return {'experiment': exp,
            'pageTitle': strings.upload_annotations_page_title,
            'formrenderer': forms.FormRenderer(schema),
            'errors':[]}