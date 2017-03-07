from pyramid.view import view_config
from ptmscout.database import experiment, upload, jobs
from ptmscout.config import strings
from ptmscout.utils import webutils, decorators, wizard
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


def create_nav_wizard(request, exp, session):
    navigation = wizard.WizardNavigation(request)

    navigation.add_page('configure_annotations', "Configure Annotation File", True, id=exp.id, sid=session.id)
    navigation.add_page('confirm_annotations', "Confirm Upload", True, id=exp.id, sid=session.id)
    navigation.set_page('confirm_annotations')

    return navigation


def upload_confirm(request, exp, session):
    confirm = webutils.post(request, "confirm", "false") == "true"
    terms_of_use_accepted = 'terms_of_use' in request.POST
    
    if session.stage == 'complete':
        raise AnnotationUploadAlreadyStarted()
    
    reason = None
    
    if confirm and terms_of_use_accepted:
        create_dataset_job(session, request)
        
    
        return {'pageTitle': strings.annotation_upload_started_page_title,
                'message': strings.annotation_upload_started_message % (request.route_url('my_experiments')),
                'experiment': exp,
                'session_id':session.id,
                'reason':reason,
                'confirm':confirm}
    
    if confirm and not terms_of_use_accepted:
        reason = strings.failure_reason_terms_of_use_not_accepted
        confirm = False
    
    return {'pageTitle': strings.annotation_upload_confirm_page_title,
            'message': strings.annotation_upload_confirm_message,
            'navigation':create_nav_wizard(request, exp, session),
            'experiment': exp,
            'session_id': session.id,
            'reason':reason,
            'confirm': confirm}
    
@view_config(route_name='confirm_annotations', renderer='ptmscout:/templates/upload/upload_confirm.pt', permission='private')
@decorators.get_experiment('id', types=set(['experiment', 'dataset']))
@decorators.get_session('sid', 'annotations')
def upload_confirm_view(context, request, exp, session):
    return upload_confirm(request, exp, session)
