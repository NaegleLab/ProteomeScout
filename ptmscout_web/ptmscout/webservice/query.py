from pyramid.view import view_config
from ptmworker import entrez_query

@view_config(route_name='pmid_fetch', renderer='json', permission='private')
def call_get_pubmed_record_by_id(request):
    pmid = int(request.matchdict['id'])
    record = entrez_query.get_pubmed_record_by_id(pmid)
    return record