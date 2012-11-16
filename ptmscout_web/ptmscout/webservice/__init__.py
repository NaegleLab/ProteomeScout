
def add_views(config):
    config.add_route('pmid_fetch', '/webservice/pubmed/{id}')
    config.add_route('field_fetch', '/webservice/autocomplete/{field}')