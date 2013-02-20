from pyramid.view import view_config
import ptmscout.database.experiment as experiment
from ptmscout.config import strings


@view_config(route_name='experiments', renderer='ptmscout:templates/experiments/experiments.pt')
def experiment_listing(request):
    experiments = experiment.getExperimentTree(request.user)
    
    return {'pageTitle': strings.experiments_page_title,
            'experiments': experiments}

@view_config(route_name='experiment', renderer='ptmscout:templates/experiments/experiment_home.pt')
def view_experiment(request):
    experiment_id = request.matchdict['id']
    ptm_exp = experiment.getExperimentById(experiment_id, request.user)
        
    return {'pageTitle': strings.experiment_page_title % (ptm_exp.name),
            'experiment': ptm_exp}

