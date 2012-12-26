from ptmscout.config import strings
from ptmscout.database import protein, modifications
from pyramid.view import view_config


@view_config(route_name='protein_mod_sites', renderer='ptmscout:templates/proteins/protein_modifications.pt')
def protein_modifications_view(request):
    pid = int(request.matchdict['id'])
    prot = protein.getProteinById(pid)
    
    mod_sites = {}
    
    for MS in modifications.getMeasuredPeptidesByProtein(pid, request.user):
        for pepmod in MS.peptides:
            pep = pepmod.peptide
            pep_tuple = (pep, pepmod.modification.name)
            exps = mod_sites.get(pep_tuple, set())
            exps.add(MS.experiment)
            mod_sites[pep_tuple] = exps
            
    mod_sites = [{'site':pep.site_pos, 'name':pep.getName(), 'type':mod_type, 'peptide':pep.getPeptide(), 'experiments':exps} for (pep, mod_type), exps in mod_sites.items()]
    mod_sites = sorted( mod_sites, key=lambda item: item['site'] )
    
    return {'protein': prot,
            'pageTitle': strings.protein_modification_sites_page_title,
            'modification_sites': mod_sites}
