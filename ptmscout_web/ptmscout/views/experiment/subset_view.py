from pyramid.view import view_config
from ptmscout.database import experiment
from ptmscout.config import strings
from ptmscout.views.dataset import dataset_explorer_view

@view_config(route_name='experiment_subset', renderer='ptmscout:templates/experiments/experiment_subset.pt')
def view_experiment_subset(request):
    experiment_id = request.matchdict['id']
    ptm_exp = experiment.getExperimentById(experiment_id, request.user)

    
    result = { 'pageTitle': strings.experiment_subset_page_title,
               'experiment': ptm_exp }
    result.update( dataset_explorer_view.format_explorer_view( ptm_exp.measurements ) )
    
    return result
