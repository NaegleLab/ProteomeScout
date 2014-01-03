from pyramid.view import view_config
from ptmscout.database import experiment
from ptmscout.config import strings
from ptmscout.views.dataset import dataset_explorer_view
from ptmscout.utils import webutils

@view_config(route_name='experiment_subset', renderer='ptmscout:templates/experiments/experiment_subset.pt', permission='private')
def view_experiment_subset(request):
    experiment_id = int(request.matchdict['id'])
    annotation_set_id = webutils.get_field_as_int(request, 'annotation_set', None)

    ptm_exp = experiment.getExperimentById(experiment_id, request.user)

    initial_query = webutils.get(request, 'query', "")
    
    result = { 'pageTitle': strings.experiment_subset_page_title % (ptm_exp.name),
               'experiment': ptm_exp,
               'initial_query': initial_query }
    result.update( dataset_explorer_view.format_explorer_view( ptm_exp, annotation_set_id, request.user ) )
    
    return result
