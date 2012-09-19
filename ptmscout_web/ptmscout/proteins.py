from ptmscout import strings
from pyramid.view import view_config
from ptmscout.database import protein, modifications
from ptmscout.utils import webutils


@view_config(route_name='protein_mod_sites', renderer='templates/protein_modifications.pt')
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
    
@view_config(route_name='protein_main', renderer='templates/protein_modifications.pt')
def protein_view(request):
    return protein_modifications_view(request)

@view_config(route_name='protein_search', renderer='templates/protein_search.pt')
def protein_search_view(request):
    submitted = bool(webutils.post(request, 'submitted', False))
    acc_search = webutils.post(request, 'acc_search', "")
    stringency = webutils.post(request, 'stringency', "1")
    selected_species = webutils.post(request, 'species', None)
    
    proteins = []
    protein_mods={}
    if(submitted and acc_search != ""):
        proteins = protein.getProteinsByAccession([acc_search], species=selected_species)
        
        for p in proteins:
            mods = modifications.getModificationsByProtein(p.id, request.user)
            for mod in mods:
                peps = protein_mods.get(p.id, set())
                [ peps.add(pep) for pep in mod.phosphopeps ]
                protein_mods[p.id] = peps
    
    proteins = sorted(proteins, key=lambda prot: prot.acc_gene)
    
    for pid in protein_mods:
        protein_mods[pid] = sorted(protein_mods[pid], key=lambda mod: mod.getName())
        protein_mods[pid] = [ { 'site': pep.getName(), 'peptide': pep.getPeptide() } for pep in protein_mods[pid]]
    
    species_list = [ species.name for species in protein.getAllSpecies() ]
    
    return {'pageTitle': strings.protein_search_page_title,
            'species_list': species_list,
            'acc_search':acc_search,
            'stringency':stringency,
            'selected_species':selected_species,
            'proteins':proteins,
            'modifications':protein_mods,
            'submitted': submitted}