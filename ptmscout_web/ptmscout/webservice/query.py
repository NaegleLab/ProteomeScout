from pyramid.view import view_config
from ptmworker.helpers import entrez_tools
from ptmscout.database import experiment
from ptmscout.utils import webutils

@view_config(route_name='pmid_fetch', renderer='json', permission='private')
def call_get_pubmed_record_by_id(request):
    pmid = int(request.matchdict['id'])
    record = entrez_tools.get_pubmed_record_by_id(pmid)
    return record

@view_config(route_name='field_fetch', renderer='json', permission='private')
def get_autocomplete_for_field(request):
    field_name = request.matchdict['field']
    field_values = experiment.getValuesForField(field_name)
    field_values = sorted(field_values, key=lambda item: item.lower())
    
    return {field_name:field_values}


def format_experiments(exps):
    return 

@view_config(route_name='experiment_fetch', renderer='json', permission='private')
def get_matching_experiments_for_user(request):
    search_term = webutils.get(request, 'search_term', None)
    conditions = webutils.get(request, 'conditions', None)
    offset = int( webutils.get(request, 'offset', 0) )
    limit = 10

    cond_map = {}
    try:
        for kv in conditions.split(','):
            k, v = kv.split(':')
            k, v = k.strip(), v.strip()
            vmap = cond_map.get(k, [])
            vmap.append(v)
            cond_map[k] = vmap
    except:
        cond_map = {}

    page = limit, offset

    exp_url = request.application_url + "/experiments/%d"

    count, experiments = experiment.searchExperiments(text_search=search_term, conditions = cond_map, user=request.user, page=page)
    formatted_experiments = [ { 'id':e.id, 'link':exp_url % (e.id), 'name':e.name, 'residues':e.modified_residues } for e in experiments ]

    return {'count': count, 'offset': offset, 'limit': limit, 'experiments': formatted_experiments}
