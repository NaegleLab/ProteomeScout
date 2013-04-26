from pyramid.view import view_config
from ptmscout.database import experiment, upload, modifications, jobs
from ptmscout.config import strings
from ptmscout.utils import webutils, decorators
from ptmworker import data_import

class UploadAlreadyStarted(Exception):
    pass
    
@view_config(context=UploadAlreadyStarted, renderer='ptmscout:/templates/info/information.pt')
def upload_already_started_view(request):
    return {'pageTitle': strings.experiment_upload_started_page_title,
            'header': strings.experiment_upload_started_page_title,
            'message': strings.experiment_upload_started_message % (request.application_url + "/account/experiments")}



def prepare_experiment(session, exp, user):
    exp_target = exp
    
    if session.load_type=='reload' or session.load_type == 'append':
        parent_exp = experiment.getExperimentById(session.parent_experiment, user, check_ready=False)
        parent_exp.copyData(exp)

        session.experiment_id = parent_exp.id
        exp_target = parent_exp

    if session.load_type=='reload':
        modifications.deleteExperimentData(parent_exp.id)
    
    exp_target.saveExperiment()
    
    session.stage = 'complete'
    session.save()
    
    if session.load_type=='reload' or session.load_type == 'append':
        exp.delete()
        
    return exp_target

def create_job(request, exp, session, user):
    job = jobs.Job()
    job.name = "Load experiment '%s'" % (exp.name)
    job.type = 'load_experiment'
    job.user_id = user.id
    job.status_url = request.route_url('my_experiments')
    job.resume_url = request.route_url('upload_resume', id=session.id)
    job.result_url = request.route_url('experiment', id=exp.id)
    job.status = 'in queue'
    job.stage = 'query'
    job.save()
    
    exp.job_id = job.id
    exp.saveExperiment()
    
    return job.id

def upload_confirm(request, session):
    confirm = webutils.post(request, "confirm", "false") == "true"
    terms_of_use_accepted = 'terms_of_use' in request.POST
    
    if session.stage == 'complete':
        raise UploadAlreadyStarted()
    
    reason = None
    
    exp = experiment.getExperimentById(session.experiment_id, request.user, False)
    exp_dict = webutils.object_to_dict(exp)
    exp_dict['citation'] = exp.getLongCitationString()
    exp_dict['url'] = exp.getUrl()
    
    if confirm and terms_of_use_accepted:
        target_exp = prepare_experiment(session, exp, request.user)
        job_id = create_job(request, target_exp, session, request.user)
        data_import.start_import.apply_async((target_exp.id, session.id, job_id))
        
        return {'pageTitle': strings.experiment_upload_started_page_title,
                'message': strings.experiment_upload_started_message % (request.application_url + "/account/experiments"),
                'experiment': exp_dict,
                'session_id':session.id,
                'reason':reason,
                'confirm':confirm}
    
    if confirm and not terms_of_use_accepted:
        reason = strings.failure_reason_terms_of_use_not_accepted
        confirm = False
    
    return {'pageTitle': strings.experiment_upload_confirm_page_title,
            'message': strings.experiment_upload_confirm_message,
            'experiment': exp_dict,
            'session_id': session.id,
            'reason':reason,
            'confirm': confirm}



@view_config(route_name='upload_confirm', renderer='ptmscout:/templates/upload/upload_confirm.pt', permission='private')
@decorators.get_session('id', 'experiment')
def upload_confirm_view(context, request, session):
    return upload_confirm(request, session)
