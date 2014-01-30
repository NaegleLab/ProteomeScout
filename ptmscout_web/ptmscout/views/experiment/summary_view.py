from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.database import experiment, modifications
import base64
import json
from ptmscout.utils import protein_utils, decorators

def summarize_measurements(measurements):
    summary = {'modifications':0,
               'measured':0,
               'proteins':set(),
               'by_residue':{},
               'by_species':{},
               'by_type':{}}
    
    for measured_peptide in measurements:
        for p in measured_peptide.peptides:
            pep = p.peptide
            
            residue_count = summary['by_residue'].get(pep.site_type, 0)
            summary['by_residue'][pep.site_type] = residue_count+1
            
            species = measured_peptide.protein.species.name
            species_count = summary['by_species'].get(species, 0)
            summary['by_species'][species] = species_count+1
            
            mod = p.modification
            while mod.parent:
                mod = mod.parent 
            
            type_count = summary['by_type'].get(mod.name, 0)
            summary['by_type'][mod.name] = type_count+1
            
            mods = summary['modifications']
            summary['modifications'] = mods+1
            
        measured = summary['measured']
        summary['measured'] = measured+1
        
        summary['proteins'].add(measured_peptide.protein_id)
    
    summary['proteins'] = len(summary['proteins'])
    
    return summary

@decorators.cache_result
def summarize_experiment(exp):
    measurement_summary = summarize_measurements(exp.measurements)
    sequence_profile = protein_utils.create_sequence_profile(exp.measurements)
    rejected_peps = len(set([err.peptide for err in exp.errors]))

    return measurement_summary, sequence_profile, rejected_peps

@view_config(route_name='experiment_summary', renderer='ptmscout:templates/experiments/experiment_summary.pt')
@decorators.get_experiment('id',types=set(['experiment','dataset']))
def experiment_summary_view(context, request, exp):
    user_owner = request.user != None and request.user.experimentOwner(exp)

    result = summarize_experiment(exp)
    print "Tuple size: %d" % (len(result))

    measurement_summary, sequence_profile, rejected_peps = result
    
    encoded = base64.b64encode(json.dumps(sequence_profile))
    return {'experiment':exp,
            'measurement_summary':measurement_summary,
            'sequence_profile': encoded,
            'rejected_peptides':rejected_peps,
            'pageTitle': strings.experiment_summary_page_title % (exp.name),
            'user_owner': user_owner}
