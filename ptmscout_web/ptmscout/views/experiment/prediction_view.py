from pyramid.view import view_config
from ptmscout.database import experiment, modifications
from ptmscout.config import strings
import base64
import json
from ptmscout.utils import decorators

def filter_predictions(predictions, threshold = 1.0):
    return [ p for p in predictions if p.percentile <= threshold ]

def format_predictions(measurements):
    formatted_predictions = {}
    
    predictions = []
    phosphomap = {}
    for m in measurements:
        for p in m.peptides:
            pep = p.peptide
            
            for s in pep.predictions:
                predictions.append(s)
                phosphomap[s] = m
    
    predictions = filter_predictions(predictions)
    
    for scansite in predictions:
        formatted_predictions[scansite.source] = {}
        
    for scansite in predictions:
        m = phosphomap[scansite]
        
        measureset = formatted_predictions[scansite.source].get(scansite.value, set())
        measureset.add(m)
        formatted_predictions[scansite.source][scansite.value] = measureset
        
    allmeasures = set(measurements)
    keyset = formatted_predictions.keys()[:]
    
    for source in keyset:
        
        total = set()
        for value in formatted_predictions[source]:
            [ total.add(m) for m in formatted_predictions[source][value] ]
            formatted_predictions[source][value] = len(formatted_predictions[source][value])
            
        diffset = allmeasures - total
        
        formatted_predictions[source]["None"] = len(diffset)
        
        formatted_predictions[source] = sorted(formatted_predictions[source].items(), key=lambda item: -item[1])
        table = formatted_predictions[source]
        
        jsontable = base64.b64encode(json.dumps(table))
        
        formatted_predictions[source] = {'json': jsontable, 'table':table}
        
        if source in strings.prediction_type_map:
            tmp = formatted_predictions[source]
            del formatted_predictions[source]
            formatted_predictions[strings.prediction_type_map[source]] = tmp 
        
    return formatted_predictions

def create_query_generator(field):
    from ptmscout.utils.query_generator import generate_scansite_query
    def query_generator(value):
        return {'query': generate_scansite_query(field, value)}

    return query_generator

@decorators.cache_result
def wrap_format_predictions(exp):
    return format_predictions(exp.measurements)

@view_config(route_name='experiment_predictions', renderer='ptmscout:templates/experiments/experiment_predictions.pt')
@decorators.get_experiment('id',types=set(['experiment']))
def prediction_view(context, request, exp):
    user_owner = request.user != None and request.user.experimentOwner(exp)

    formatted_predictions = wrap_format_predictions(exp)

    return {'experiment': exp,
            'pageTitle': strings.experiment_prediction_page_title % (exp.name),
            'predictions': formatted_predictions,
            'user_owner': user_owner,
            'query_generators': {'Scansite Bind': create_query_generator('Scansite-Bind'), 'Scansite Kinase': create_query_generator('Scansite-Kinase')}
            }
