from pyramid.view import view_config
from ptmscout.database import experiment
from pyramid.httpexceptions import HTTPFound
from ptmscout.config import strings
from ptmscout.utils import webutils
from ptmworker import tasks

class UploadAlreadyStarted(Exception):
    def __init__(self,eid):
        self.eid = eid
    
@view_config(context=UploadAlreadyStarted, renderer='ptmscout:/templates/info/information.pt')
def upload_already_started_view(request):
    return {'pageTitle': strings.experiment_upload_started_page_title,
            'header': strings.experiment_upload_started_page_title,
            'message': strings.experiment_upload_started_message % (request.application_url + "/account/experiments")}

@view_config(route_name='upload_status', renderer='ptmscout:/templates/experiments/upload_confirm.pt')
def upload_confirm_view(request):
    confirm = webutils.post(request, "confirm", "false")
    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(int(eid), request.user, False)
    
    if exp.ready():
        return HTTPFound(request.application_url + "/experiment/" + str(eid))

    if exp.status == 'loading':
        raise UploadAlreadyStarted(eid)
    
    if confirm == "true":
        tasks.start_import.apply_async(exp)
        return {'pageTitle': strings.experiment_upload_started_page_title,
                'message': strings.experiment_upload_started_message % (request.application_url + "/account/experiments"),
                'experiment': exp,
                'confirm':confirm}
    
    return {'pageTitle': strings.experiment_upload_confirm_page_title,
            'message': strings.experiment_upload_confirm_message,
            'experiment': exp,
            'confirm': confirm}