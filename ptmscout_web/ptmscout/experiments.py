from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound, HTTPForbidden
import database.experiment as experiment
from ptmscout import strings

@view_config(route_name='redirect_to_experiments')
def redirect_to_experiments(request):
    return HTTPFound(request.application_url + "/experiments")

@view_config(route_name='experiments', renderer='templates/experiments.pt')
def experiment_listing(request):
    experiments = experiment.getExperimentTree()
    
    permitted_ids = []
    if request.user != None:
        permitted_ids = [p.experiment_id for p in request.user.permissions]
        
    experiments = [exp for exp in experiments if (exp.public == 1 or exp.id in permitted_ids)]
    
    return {'pageTitle': strings.experiments_page_title,
            'experiments': experiments}

@view_config(route_name='experiment', renderer='templates/experiment_home.pt')
def view_experiment(request):
    experiment_id = request.matchdict['id']
    
    ptm_exp = experiment.getExperimentById(experiment_id)
    
    if not ptm_exp.public and request.user == None:
        raise HTTPForbidden()
    
    if not ptm_exp.public and request.user != None:
        permitted_ids = [p.experiment_id for p in request.user.permissions]
        if ptm_exp.id not in permitted_ids:
            raise HTTPForbidden() 
        
    return {'pageTitle': strings.experiment_page_title % (ptm_exp.name),
            'experiment': ptm_exp}