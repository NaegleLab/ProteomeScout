from ptmscout.database import upload, experiment, modifications
from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound
from ptmscout.views.upload import upload_confirm
from ptmscout.utils import webutils, decorators
from ptmworker import data_import
from ptmscout.config import strings


def resume_upload_session(request, session):
    base_url = request.application_url + "/upload/%d/%s"
    
    if session.stage == 'complete':
        return HTTPFound(base_url % (session.id, 'confirm'))
    elif session.stage == 'condition':
        return HTTPFound(base_url % (session.id, 'conditions'))
    else:
        return HTTPFound(base_url % (session.id, session.stage))

@view_config(route_name='upload_resume', permission='private')
@decorators.get_session('id', 'experiment')
def resume_upload_session_view(context, request, session):
    return resume_upload_session(request, session)

def prepare_experiment(session, exp):
    exp.job.status='in queue'
    exp.job.restart()
    exp.job.save()
    exp.saveExperiment()

def retry_failed_upload(request, session):
    exp = experiment.getExperimentById(session.experiment_id, request.user, False)

    if exp.status != 'error':
        raise upload_confirm.UploadAlreadyStarted()
    
    prepare_experiment(session, exp)
    data_import.start_import.apply_async((exp.id, session.id, exp.job.id))

    return {'pageTitle': strings.experiment_upload_started_page_title,
            'message': strings.experiment_upload_started_message % (request.application_url + "/account/experiments"),
            'experiment': exp,
            'session_id':session.id,
            'reason':None,
            'confirm':True}

@view_config(route_name='upload_retry', renderer='ptmscout:/templates/upload/upload_confirm.pt', permission='private')
@decorators.get_session('id', 'experiment')
def retry_failed_upload_view(context, request, session):
    return retry_failed_upload(request, session)
