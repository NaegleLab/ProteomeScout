from ptmscout.config import strings
from ptmscout.database import protein
from pyramid.view import view_config
from ptmscout.views.protein import decorators


@view_config(route_name='protein_GO', renderer='ptmscout:templates/proteins/protein_ontology.pt')
@decorators.experiment_filter
def protein_gene_ontology_view(context, request):
    pid = int(request.matchdict['id'])
    prot = protein.getProteinById(pid)
    
    term_dict = {'F':set(),'P':set(),'C':set()}
    
    for goe in prot.GO_terms:
        term = goe.GO_term
        term_dict[term.aspect].add((term.GO, term.term))
        
    for aspect in term_dict:
        term_dict[aspect] = sorted(list(term_dict[aspect]), key=lambda item: item[1])
    
    return {'pageTitle': strings.protein_ontology_page_title,
            'protein': prot,
            'F_terms':term_dict['F'],
            'P_terms':term_dict['P'],
            'C_terms':term_dict['C']}
