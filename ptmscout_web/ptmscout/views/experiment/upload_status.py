from pyramid.view import view_config
from ptmscout.database import experiment
from pyramid.httpexceptions import HTTPFound
from ptmscout.config import strings
from ptmscout.utils import webutils



@view_config(route_name='upload_status', renderer='ptmscout:/templates/experiments/upload_confirm.pt')
def upload_confirm_view(request):
    confirm = webutils.post(request, "confirm", "false")
    eid = request.matchdict['id']
    exp = experiment.getExperimentById(int(eid), request.user, False)
    
    if exp.ready == 1:
        return HTTPFound(request.application_url + "/experiment/" + eid)
    
    return {'pageTitle': strings.experiment_upload_confirm_page_title,
            'message': strings.experiment_upload_confirm_message,
            'experiment': exp,
            'confirm': confirm}