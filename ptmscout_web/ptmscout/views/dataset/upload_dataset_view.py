from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound
from ptmscout.config import settings, strings
from ptmscout.utils import uploadutils, forms
from ptmscout.database import experiment, upload

def create_experiment(name, filename, user):
    exp = experiment.Experiment()
    exp.name = name
    exp.description=''
    exp.published=0
    exp.ambiguity=0
    exp.public=0
    exp.dataset=filename
    exp.type='dataset'

    exp.grantPermission(user, 'owner')
    exp.saveExperiment()

    return exp.id

def create_session(experiment_id, annotation_filename, user):
    session = upload.Session()
    session.data_file = annotation_filename
    session.experiment_id = experiment_id
    
    session.load_type = 'dataset'
    session.stage = 'config'
    session.user_id = user.id
    session.change_name=''
    session.change_description = ""
    
    session.save()

    return session.id

def create_schema(request):
    schema = forms.FormSchema()
    
    schema.add_text_field('datasetname', 'Dataset Name')
    schema.add_file_upload_field('datasetfile', 'Input Data File')
    schema.set_required_field('datasetname')
    schema.set_required_field('datasetfile')
    schema.parse_fields(request)
    
    return schema

@view_config(route_name='dataset_upload', request_method='POST', renderer='ptmscout:/templates/dataset/upload_dataset.pt', permission='private')
def upload_dataset_POST(request):
    schema = create_schema(request)
    errors = forms.FormValidator(schema).validate()
    
    if len(errors) == 0:
        output_file = uploadutils.save_data_file(request.POST['datasetfile'], settings.dataset_files_prefix)
        experiment_id = create_experiment(request.POST['datasetname'], output_file, request.user)
        sid = create_session(experiment_id, output_file, request.user)

        return HTTPFound(request.route_url('dataset_configure', id=sid))

    return {
            'pageTitle': strings.upload_annotations_page_title,
            'formrenderer': forms.FormRenderer(schema),
            'errors':errors}
    
@view_config(route_name='dataset_upload', request_method='GET', renderer='ptmscout:/templates/dataset/upload_dataset.pt', permission='private')
def upload_dataset_GET(request):
    schema = create_schema(request)
    
    return {
            'pageTitle': strings.upload_annotations_page_title,
            'formrenderer': forms.FormRenderer(schema),
            'errors':[]}
