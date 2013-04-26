from pyramid.view import view_config
from ptmscout.database import experiment, jobs
from ptmscout.config import strings
from ptmscout.utils import webutils, decorators
from ptmworker import data_import

class DatasetUploadAlreadyStarted(Exception):
    pass
    
@view_config(context=DatasetUploadAlreadyStarted, renderer='ptmscout:/templates/info/information.pt')
def upload_already_started_view(request):
    return {'pageTitle': strings.dataset_upload_started_page_title,
            'header': strings.dataset_upload_started_page_title,
            'message': strings.dataset_upload_started_message % (request.application_url + "/account/experiments")}


def create_dataset_job(exp, session, request):
    job = jobs.Job()
    job.name = 'Load Dataset: %s' % (exp.name)
    job.type = 'load_dataset'
    
    job.stage = 'query'
    job.status = 'in queue'
    job.status_url = request.route_url('my_experiments')
    job.result_url = request.route_url('experiment', id=exp.id)
    job.user_id = request.user.id
    
    job.save()
    
    exp.job_id = job.id
    exp.saveExperiment()
        
    data_import.start_import.apply_async((exp.id, session.id, job.id, True))

def upload_confirm(request, session):
    exp = experiment.getExperimentById(session.experiment_id, request.user, False)
    exp_dict = webutils.object_to_dict(exp)
    exp_dict['citation'] = exp.getLongCitationString()
    exp_dict['url'] = exp.getUrl()
    
    confirm = webutils.post(request, "confirm", "false") == "true"
    terms_of_use_accepted = 'terms_of_use' in request.POST
    
    if session.stage == 'complete':
        raise DatasetUploadAlreadyStarted()
    
    reason = None
    
    if confirm and terms_of_use_accepted:
        create_dataset_job(exp, session, request)
        
    
        return {'pageTitle': strings.dataset_upload_started_page_title,
                'message': strings.dataset_upload_started_message % (request.route_url('my_experiments')),
                'experiment': exp_dict,
                'session_id':session.id,
                'reason':reason,
                'confirm':confirm}
    
    if confirm and not terms_of_use_accepted:
        reason = strings.failure_reason_terms_of_use_not_accepted
        confirm = False
    
    return {'pageTitle': strings.dataset_upload_confirm_page_title,
            'message': strings.dataset_upload_confirm_message,
            'experiment': exp_dict,
            'session_id': session.id,
            'reason':reason,
            'confirm': confirm}
    
    
@view_config(route_name='dataset_confirm', renderer='ptmscout:/templates/upload/upload_confirm.pt', permission='private')
@decorators.get_session('id', 'dataset')
def upload_confirm_view(context, request, session):
    return upload_confirm(request, session)