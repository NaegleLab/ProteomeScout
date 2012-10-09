from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.database import experiment, modifications
import base64
import json

def build_go_annotation_tree(measurements):
    tree = {'F':[], 'P':[], 'C':[], 'total':0}
    
    protein_set = set()
    for m in measurements:
        protein_set.add(m.protein)
    
    GO_terms = {}
    for p in protein_set:
        for g in p.GO_terms:
            node = GO_terms.get(g.GO, {'GO':g.GO, 'aspect':g.aspect, 'term':g.term, 'value':0, 'children':[]})
            node['value'] += 1
            GO_terms[g.GO] = node
            
    GO_set = set(GO_terms.keys())
    child_set = set()
    
    for p in protein_set:
        for g in p.GO_terms:
            parent = GO_terms[g.GO]
            for c in g.children:
                if c.GO in GO_terms: 
                    child_set.add(c.GO)
                    child = GO_terms[c.GO]
                    if(child not in parent['children']):
                        parent['children'].append(child)
    
    for GO in GO_terms:
        GO_terms[GO]['children'] = sorted(GO_terms[GO]['children'], key=lambda item: item['GO'])
    
    root_set = GO_set - child_set
    GO_roots = []    
    
    for GO in root_set:
        GO_roots.append(GO_terms[GO])
        
    GO_roots = sorted(GO_roots, key=lambda item: item['GO'])
    
    tree['F'] = [ term for term in GO_roots if term['aspect'] == 'F' ]
    tree['P'] = [ term for term in GO_roots if term['aspect'] == 'P' ]
    tree['C'] = [ term for term in GO_roots if term['aspect'] == 'C' ]
    
    tree['total'] = sum([ term['value'] for term in GO_roots ])
    
    return base64.b64encode(json.dumps(tree))
    

def format_go_terms(measurements):
    GO_terms = {'F':{},'P':{},'C':{}}
    prot_by_aspect = {'F':set(), 'P':set(), 'C':set()}
    
    protein_set = set()
    for m in measurements:
        protein_set.add(m.protein)
    
    for p in protein_set:
        for g in p.GO_terms:
            num = GO_terms[g.aspect].get((g.GO, g.term), 0)
            GO_terms[g.aspect][(g.GO, g.term)] = num+1
            prot_by_aspect[g.aspect].add(p)
            
    def comparator(((G1, T1), V1), ((G2, T2), V2)):
        if V1 == V2:
            return -1 if G1 < G2 else (0 if G1 == G2 else 1)
        return V2 - V1
    
    for aspect in GO_terms:
        GO_terms[aspect][('None', '-')] = len(protein_set - prot_by_aspect[aspect])
        
        terms = GO_terms[aspect].items()
        terms.sort(comparator)
        GO_terms[aspect] = [(GO, term, cnt) for ((GO, term), cnt) in terms]
    
    
    return {'molecular_function':GO_terms['F'],
              'cellular_component':GO_terms['C'],
              'biological_process':GO_terms['P']}

@view_config(route_name='experiment_GO', renderer='ptmscout:/templates/experiments/experiment_GO.pt')
def show_experiment_go_terms(request):
    exp_id = int(request.matchdict['id'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    
    measurements = modifications.getMeasuredPeptidesByExperiment(exp_id, request.user)
    formatted_go_terms = format_go_terms(measurements)
    go_tree = build_go_annotation_tree(measurements)
    
    return {'pageTitle': strings.experiment_GO_page_title % (exp.name),
            'experiment': exp,
            'go_tables': formatted_go_terms,
            'go_tree': go_tree}