from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound, HTTPForbidden
import database.experiment as experiment
from ptmscout import strings

@view_config(route_name='redirect_to_experiments')
def redirect_to_experiments(request):
    return HTTPFound(request.application_url + "/experiments")

@view_config(route_name='experiments', renderer='templates/experiments.pt')
def experiment_listing(request):
    experiments = experiment.getExperimentTree(request.user)
    
    return {'pageTitle': strings.experiments_page_title,
            'experiments': experiments}

@view_config(route_name='experiment', renderer='templates/experiment_home.pt')
def view_experiment(request):
    experiment_id = request.matchdict['id']
    ptm_exp = experiment.getExperimentById(experiment_id, request.user)
        
    return {'pageTitle': strings.experiment_page_title % (ptm_exp.name),
            'experiment': ptm_exp}