from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.database import experiment, modifications
import base64
import json
import math


def create_sequence_profile(measurements):
    peptides = [p.peptide for m in measurements for p in m.peptides]

    frequencies = [0]*15
    N = float(len(peptides))
    
    for i in xrange(0, 15):
        frequencies[i] = {}
       
    for pep in peptides:
        sequence = pep.pep_aligned.upper()
        
        for i, s in enumerate(sequence):
            if s == ' ':
                s = '-'
            val = frequencies[i].get(s, 0)
            frequencies[i][s] = val+1
        
    seqlogo = {'total':len(peptides), 'frequencies':[]}
    
    if len(peptides) == 0:
        return seqlogo
    
    en = 19 / (2 * math.log(2) * len(peptides)) 
    
    for i in xrange(0, 15):
        Ri = math.log(20, 2)
        
        for s in frequencies[i]:
            f = frequencies[i][s] / N
            Ri += f * math.log(f, 2)
        
        Ri -= en
        
        sorted_freqs = sorted([ (k, v) for (k,v) in frequencies[i].items() ], key=lambda item: -item[1])
        final = []
        
        d=0
        for k,v in sorted_freqs:
            final.append((k, v, d))
            d += v
                
        seqlogo['frequencies'].append({'R':Ri, 'f':final})
    
    return seqlogo

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

@view_config(route_name='experiment_summary', renderer='ptmscout:templates/experiments/experiment_summary.pt')
def experiment_summary_view(request):
    eid = request.matchdict['id']
    exp = experiment.getExperimentById(eid, request.user)
    measurements = modifications.getMeasuredPeptidesByExperiment(eid, request.user)
    
    user_owner = request.user != None and request.user.experimentOwner(exp)
    
    rejected_peps = len(set([err.peptide for err in exp.errors]))
    
    measurement_summary = summarize_measurements(measurements)
    sequence_profile = create_sequence_profile(measurements)
    
    encoded = base64.b64encode(json.dumps(sequence_profile))
    return {'experiment':exp,
            'measurement_summary':measurement_summary,
            'sequence_profile': encoded,
            'rejected_peptides':rejected_peps,
            'pageTitle': strings.experiment_summary_page_title % (exp.name),
            'user_owner': user_owner}
