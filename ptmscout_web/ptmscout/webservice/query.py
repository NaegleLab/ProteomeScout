from pyramid.view import view_config
from ptmworker import entrez_tools
from ptmscout.database import experiment

@view_config(route_name='pmid_fetch', renderer='json', permission='private')
def call_get_pubmed_record_by_id(request):
    pmid = int(request.matchdict['id'])
    record = entrez_tools.get_pubmed_record_by_id(pmid)
    return record

@view_config(route_name='field_fetch', renderer='json', permission='private')
def get_autocomplete_for_field(request):
    field_name = request.matchdict['field']
    field_values = experiment.getValuesForField(field_name)
    
    return {field_name:field_values}