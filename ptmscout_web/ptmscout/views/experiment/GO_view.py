from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.database import experiment, modifications

def build_go_annotation_tree(measurements):
    return ""
    

def format_go_terms(measurements):
    GO_terms = {'F':{},'P':{},'C':{}}
    
    protein_set = set()
    for m in measurements:
        protein_set.add(m.protein)
            
    for p in protein_set:
        for g in p.GO_terms:
            num = GO_terms[g.aspect].get((g.GO, g.term), 0)
            GO_terms[g.aspect][(g.GO, g.term)] = num+1
    
    for aspect in GO_terms:
        sorted_terms = sorted(GO_terms[aspect].items(), key=lambda item:-item[1])
        GO_terms[aspect] = [(GO, term, cnt) for ((GO, term), cnt) in sorted_terms]
    
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