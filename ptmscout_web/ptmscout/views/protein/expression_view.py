from ptmscout.config import strings
from ptmscout.database import protein
from pyramid.view import view_config


@view_config(route_name='protein_expression', renderer='ptmscout:templates/proteins/protein_expression.pt')
def protein_expression_view(request):
    pid = int(request.matchdict['id'])
    prot = protein.getProteinById(pid)
    
    probe_ids = []
    collections = set()
    expression_data = []
    
    for probe in prot.expression_probes:
        probe_ids.append(probe.probeset_id)
        
        for sample in probe.samples:
            col_name = sample.collection.name
            tissue_name = sample.tissue.name
            collections.add(col_name)
            expression_data.append({'probeset':probe.probeset_id, 
                                    'collection':col_name, 
                                    'tissue': tissue_name,
                                    'value': sample.value})

    expression_data = sorted(expression_data, key=lambda item: item['tissue'])

    return {'pageTitle': strings.protein_expression_page_title,
            'protein': prot,
            'probe_ids': sorted(probe_ids),
            'collections': sorted(list(collections)),
            'expression_data':expression_data}
