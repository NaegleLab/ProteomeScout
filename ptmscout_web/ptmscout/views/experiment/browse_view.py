from ptmscout.config import strings
from ptmscout.database import experiment, protein, modifications
from ptmscout.utils import webutils
from pyramid.view import view_config
@view_config(route_name='experiment_browse', renderer='ptmscout:templates/experiments/experiment_browse.pt')
def browse_experiment(request):
    submitted = webutils.post(request, 'submitted', False)
    acc_search = webutils.post(request, 'acc_search', "").strip()
    stringency = webutils.post(request, 'stringency', "1").strip()
    submitted = (submitted == "true")
    
    experiment_id = request.matchdict['id']
    ptm_exp = experiment.getExperimentById(experiment_id, request.user)
    
    proteins = []
    mods = {}
    predictions = {}
    
    if(not submitted or len(acc_search) > 0):
        
        if(submitted):
            protein_list = protein.getProteinsByAccession([acc_search])
            mod_list = modifications.getMeasuredPeptidesByExperiment(ptm_exp.id, request.user, [p.id for p in protein_list])
        else:
            mod_list = modifications.getMeasuredPeptidesByExperiment(ptm_exp.id, request.user)
        
        prots = {}
        
        for mod in mod_list:
            prots[mod.protein_id] = mod.protein
            
            pep_list = mods.get(mod.protein_id, set())
            
            for p in mod.peptides:
                pep = p.peptide
                pep_list.add(pep)
            
            mods[mod.protein_id] = pep_list
        
        proteins = sorted( [ prots[pid] for pid in prots ], key=lambda prot: prot.acc_gene )
        
        for pid in mods:
            predictions[pid] = []
            if len(mods[pid])==1:
                for pep in mods[pid]:
                    predictions[pid].extend(pep.predictions)
            
            mods[pid] = [ {'site':pep.getName(), 'peptide':pep.getPeptide()} for pep in mods[pid] ]
            mods[pid] = sorted( mods[pid], key=lambda pep: pep['site'] )
            
        for pid in predictions:
            pid_predictions = predictions[pid]
            predictions[pid] = {}
            
            for scansite in pid_predictions:
                predictions[pid][scansite.source] = []
                
            for scansite in pid_predictions:
                predictions[pid][scansite.source].append((scansite.value, scansite.score)) 
                
            for source in predictions[pid]:
                predictions[pid][source] = sorted(predictions[pid][source], key=lambda item: item[0])           
            
    
    return {'acc_search': acc_search,
            'stringency': stringency,
            'experiment': ptm_exp,
            'submitted': submitted,
            'pageTitle': strings.experiment_browse_page_title % (ptm_exp.name),
            'proteins': proteins,
            'include_predictions': True,
            'scansites': predictions,
            'modifications': mods}