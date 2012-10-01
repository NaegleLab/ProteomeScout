from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.database import experiment


@view_config(route_name='experiment_summary', renderer='ptmscout:/templates/experiments/experiment_summary.pt')
def experiment_summary_view(request):
    eid = request.matchdict['id']
    exp = experiment.getExperimentById(eid, request.user)
    
    return {'experiment':exp,
            'pageTitle': strings.experiment_summary_page_title % (exp.name)}