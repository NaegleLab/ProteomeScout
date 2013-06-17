from ptmscout.config import strings
from ptmscout.database import protein, modifications
from pyramid.view import view_config
from ptmscout.views.protein import decorators

@view_config(route_name='protein_mod_sites', renderer='ptmscout:templates/proteins/protein_modifications.pt')
@decorators.experiment_filter
def protein_modifications_view(context, request):
    pid = int(request.matchdict['id'])
    prot = protein.getProteinById(pid)
    
    mod_sites = {}
    experiment_filter = request.urlfilter.get_field('experiment_id')
    mspeps = modifications.getMeasuredPeptidesByProtein(pid, request.user)

    if experiment_filter:
        mspeps = [ ms for ms in mspeps if ms.experiment_id == experiment_filter ]

    for MS in mspeps:
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
