from pyramid.view import view_config
from ptmscout.database import experiment, upload, jobs
from ptmscout.config import strings
from ptmscout.utils import webutils
from ptmworker import annotate_tasks

class AnnotationUploadAlreadyStarted(Exception):
    pass
    
@view_config(context=AnnotationUploadAlreadyStarted, renderer='ptmscout:/templates/info/information.pt')
def upload_already_started_view(request):
    return {'pageTitle': strings.annotation_upload_started_page_title,
            'header': strings.annotation_upload_started_page_title,
            'message': strings.annotation_upload_started_message % (request.application_url + "/account/experiments")}


def create_dataset_job(session, request):
    job = jobs.Job()
    job.name = 'Load Annotations for Experiment %d' % (session.experiment_id)
    job.type = 'load_annotations'
    
    job.status = 'in queue'
    job.status_url = request.route_url('my_experiments')
    job.result_url = request.route_url('experiment_subset', id=session.experiment_id)
    job.user_id = request.user.id
    
    job.save()
    
    annotate_tasks.start_annotation_import.apply_async((session.id, job.id))

@view_config(route_name='confirm_annotations', renderer='ptmscout:/templates/upload/upload_confirm.pt', permission='private')
def upload_confirm_view(request):
    session_id = int(request.matchdict['sid'])
    experiment_id = int(request.matchdict['id'])
    
    session = upload.getSessionById(session_id, request.user)
    
    exp = experiment.getExperimentById(experiment_id, request.user, False)
    exp_dict = webutils.object_to_dict(exp)
    exp_dict['citation'] = exp.getLongCitationString()
    exp_dict['url'] = exp.getUrl()
    
    confirm = webutils.post(request, "confirm", "false") == "true"
    terms_of_use_accepted = 'terms_of_use' in request.POST
    
    if session.stage == 'complete':
        raise AnnotationUploadAlreadyStarted()
    
    reason = None
    
    if confirm and terms_of_use_accepted:
        create_dataset_job(session, request)
        
    
        return {'pageTitle': strings.annotation_upload_started_page_title,
                'message': strings.annotation_upload_started_message % (request.application_url + "/account/experiments"),
                'experiment': exp_dict,
                'session_id':session_id,
                'reason':reason,
                'confirm':confirm}
    
    if confirm and not terms_of_use_accepted:
        reason = strings.failure_reason_terms_of_use_not_accepted
        confirm = False
    
    return {'pageTitle': strings.annotation_upload_confirm_page_title,
            'message': strings.annotation_upload_confirm_message,
            'experiment': exp_dict,
            'session_id': session_id,
            'reason':reason,
            'confirm': confirm}