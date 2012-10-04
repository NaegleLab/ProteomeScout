from pyramid.view import view_config
from ptmscout.database import experiment, modifications
from ptmscout.config import strings
import base64
import json

def filter_predictions(predictions, threshold = 1.0):
    return [p for p in predictions if p.score <= threshold and p.value != "~~~"]

def format_predictions(measurements):
    formatted_predictions = {}
    
    predictions = []
    phosphomap = {}
    for m in measurements:
        for p in m.phosphopeps:
            for s in p.predictions:
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

@view_config(route_name='experiment_predictions', renderer='ptmscout:templates/experiments/experiment_predictions.pt')
def prediction_view(request):
    expid = request.matchdict['id']
    exp = experiment.getExperimentById(expid, request.user)
    
    
    measurements = modifications.getMeasuredPeptidesByExperiment(expid, request.user)
    formatted_predictions = format_predictions(measurements)
    
    return {'experiment': exp,
            'pageTitle': strings.experiment_prediction_page_title % (exp.name),
            'predictions': formatted_predictions}