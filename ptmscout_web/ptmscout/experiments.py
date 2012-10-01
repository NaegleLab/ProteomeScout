from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound
import database.experiment as experiment
from ptmscout import strings
from ptmscout.utils import webutils
from ptmscout.database import modifications, protein

@view_config(route_name='redirect_to_experiments')
def redirect_to_experiments(request):
    return HTTPFound(request.application_url + "/experiments")

@view_config(route_name='experiments', renderer='templates/experiments.pt')
def experiment_listing(request):
    experiments = experiment.getExperimentTree(request.user)
    
    return {'pageTitle': strings.experiments_page_title,
            'experiments': experiments}

@view_config(route_name='experiment', renderer='templates/experiment_home.pt')
def view_experiment(request):
    experiment_id = request.matchdict['id']
    ptm_exp = experiment.getExperimentById(experiment_id, request.user)
        
    return {'pageTitle': strings.experiment_page_title % (ptm_exp.name),
            'experiment': ptm_exp}

@view_config(route_name='experiment_browse', renderer='templates/experiment_browse.pt')
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
            mod_list = modifications.getModificationsByExperiment(ptm_exp.id, request.user, [p.id for p in protein_list])
        else:
            mod_list = modifications.getModificationsByExperiment(ptm_exp.id, request.user)
        
        prots = {}
        
        for mod in mod_list:
            prots[mod.protein_id] = mod.protein
            
            pep_list = mods.get(mod.protein_id, set())
            
            for pep in mod.phosphopeps:
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