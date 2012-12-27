from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.database import upload, experiment

@view_config(route_name='upload_cancel', renderer='ptmscout:/templates/info/information.pt')
def cancel_upload_view(request):
    session_id = int(request.matchdict['id'])
    session = upload.getSessionById(session_id, request.user)

    
    if session.stage == 'complete':
        return {'pageTitle':strings.cancel_upload_successful_page_title,
                'header':strings.cancel_upload_already_started_header,
                'message':strings.cancel_upload_already_started_message}
    else:
        session.delete()
        
        if(session.experiment_id != None):
            exp = experiment.getExperimentById(session.experiment_id, request.user, check_ready=False)
            exp.delete()
        
        
        return {'pageTitle':strings.cancel_upload_successful_page_title,
                'header':strings.cancel_upload_successful_header,
                'message':strings.cancel_upload_successful_message}
