from ptmscout.database import experiment, protein, taxonomies
from ptmscout.config import strings
from ptmscout.utils import forms, paginate
from ptmscout.views.protein import search_view
from pyramid.view import view_config

def query_all(exp_id, pager):
    protein_cnt, proteins = protein.getProteinsByExperiment(exp_id, pager.get_pager_limits())
    pager.set_result_size(protein_cnt)

    return proteins

@view_config(route_name='experiment_browse', renderer='ptmscout:templates/experiments/experiment_browse.pt')
def browse_experiment(request):
    species_list = [ species.name for species in taxonomies.getAllSpecies() ]
    submitted, form_schema = search_view.build_schema(request, species_list)

    pager = paginate.Paginator(form_schema, search_view.QUERY_PAGE_LIMIT)
    pager.parse_parameters(request)

    experiment_id = int(request.matchdict['id'])
    ptm_exp = experiment.getExperimentById(experiment_id, request.user)
    
    proteins = []
    protein_metadata = {}
    errors = []

    if submitted:
        errors = search_view.build_validator(form_schema).validate()
        if len(errors) == 0:
            proteins = search_view.perform_query(form_schema, pager, experiment_id)
    else:
        proteins = query_all(experiment_id, pager)

    for p in proteins:
        search_view.get_protein_metadata(p, protein_metadata, request.user, experiment_id)

    form_renderer = forms.FormRenderer(form_schema)
    return {'pageTitle': strings.experiment_browse_page_title % (ptm_exp.name),
            'experiment': ptm_exp,
            'form':form_renderer,
            'pager': pager,
            'proteins':proteins,
            'protein_metadata':protein_metadata,
            'errors': errors,
            'submitted': submitted}
