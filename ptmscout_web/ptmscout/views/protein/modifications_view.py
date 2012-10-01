from ptmscout.config import strings
from ptmscout.database import protein, modifications
from pyramid.view import view_config


@view_config(route_name='protein_mod_sites', renderer='ptmscout:templates/proteins/protein_modifications.pt')
def protein_modifications_view(request):
    pid = int(request.matchdict['id'])
    prot = protein.getProteinById(pid)
    
    mod_sites = {}
    
    for MS in modifications.getModificationsByProtein(pid, request.user):
        for pep in MS.phosphopeps:
            pep_tuple = (pep.getName(), pep.getPeptide())
            exps = mod_sites.get(pep_tuple, set())
            exps.add(MS.experiment)
            mod_sites[pep_tuple] = exps
            
    mod_sites = [{'name':name, 'peptide':pep, 'experiments':exps} for (name, pep), exps in mod_sites.items()]
    mod_sites = sorted( mod_sites, key=lambda item: item['name'] )
    
    return {'protein': prot,
            'pageTitle': strings.protein_modification_sites_page_title,
            'modification_sites': mod_sites}