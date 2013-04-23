from ptmscout.database import upload, experiment
from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound
from ptmscout.views.upload import upload_confirm
from ptmscout.utils import webutils
from ptmworker import data_import
from ptmscout.config import strings

@view_config(route_name='dataset_upload_resume', permission='private')
def resume_upload_session(request):
    session_id = int(request.matchdict['id'])
    session = upload.getSessionById(session_id, request.user)
    
    if session.stage == 'complete':
        return HTTPFound(request.route_url('dataset_confirm', id=session_id))
    else:
        return HTTPFound(request.route_url('dataset_configure', id=session_id))

def prepare_experiment(session, exp):
    exp.job.status='in queue'
    exp.job.restart()
    exp.job.save()
    exp.saveExperiment()


@view_config(route_name='dataset_upload_retry', renderer='ptmscout:/templates/upload/upload_confirm.pt', permission='private')
def retry_failed_upload(request):
    session_id = int(request.matchdict['id'])
    session = upload.getSessionById(session_id, request.user)
    exp = experiment.getExperimentById(session.experiment_id, request.user, False)

    if exp.status != 'error':
        raise upload_confirm.UploadAlreadyStarted()
    
    exp_dict = webutils.object_to_dict(exp)
    exp_dict['citation'] = exp.getLongCitationString()
    exp_dict['url'] = exp.getUrl()
    
    prepare_experiment(session, exp)
    data_import.start_import.apply_async((exp.id, session_id, exp.job.id, True))

    return {'pageTitle': strings.dataset_upload_started_page_title,
            'message': strings.dataset_upload_started_message % (request.application_url + "/account/experiments"),
            'experiment': exp_dict,
            'session_id':session_id,
            'reason':None,
            'confirm':True}
