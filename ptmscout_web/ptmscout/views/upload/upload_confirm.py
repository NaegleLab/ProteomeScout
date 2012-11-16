from pyramid.view import view_config
from ptmscout.database import experiment, upload
from ptmscout.config import strings
from ptmscout.utils import webutils
from ptmworker import data_import

class UploadAlreadyStarted(Exception):
    pass
    
@view_config(context=UploadAlreadyStarted, renderer='ptmscout:/templates/info/information.pt')
def upload_already_started_view(request):
    return {'pageTitle': strings.experiment_upload_started_page_title,
            'header': strings.experiment_upload_started_page_title,
            'message': strings.experiment_upload_started_message % (request.application_url + "/account/experiments")}

@view_config(route_name='upload_confirm', renderer='ptmscout:/templates/upload/upload_confirm.pt')
def upload_confirm_view(request):
    confirm = webutils.post(request, "confirm", "false") == "true"
    terms_of_use_accepted = 'terms_of_use' in request.POST
    session_id = int(request.matchdict['id'])
    
    session = upload.getSessionById(session_id, request.user)
    
    if session.stage == 'complete':
        raise UploadAlreadyStarted()
    
    reason = None
    
    exp = experiment.getExperimentById(session.experiment_id, request.user, False)
    exp_dict = webutils.object_to_dict(exp)
    exp_dict['citation'] = exp.getLongCitationString()
    exp_dict['url'] = exp.getUrl()
    
    if confirm and terms_of_use_accepted:
        session.stage = 'complete'
        session.save()
        
        data_import.start_import.apply_async((exp.id, session.id, request.user.email, request.application_url))
        
        return {'pageTitle': strings.experiment_upload_started_page_title,
                'message': strings.experiment_upload_started_message % (request.application_url + "/account/experiments"),
                'experiment': exp_dict,
                'session_id':session_id,
                'reason':reason,
                'confirm':confirm}
    
    if confirm and not terms_of_use_accepted:
        reason = strings.failure_reason_terms_of_use_not_accepted
        confirm = False
    
    return {'pageTitle': strings.experiment_upload_confirm_page_title,
            'message': strings.experiment_upload_confirm_message,
            'experiment': exp_dict,
            'session_id': session_id,
            'reason':reason,
            'confirm': confirm}