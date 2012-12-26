from ptmscout.config import strings
from ptmscout.database import protein, modifications, taxonomies
from ptmscout.utils import webutils
from pyramid.view import view_config


@view_config(route_name='protein_search', renderer='ptmscout:templates/proteins/protein_search.pt')
def protein_search_view(request):
    submitted = bool(webutils.post(request, 'submitted', False))
    acc_search = webutils.post(request, 'acc_search', "")
    pep_search = webutils.post(request, 'pep_search', "")
    stringency = webutils.post(request, 'stringency', "1")
    selected_species = webutils.post(request, 'species', None)
    
    if(selected_species == 'all'):
        selected_species = None
    
    proteins = []
    protein_mods={}
    if(submitted and ( acc_search != "" or pep_search != "" )):
        protein_cnt, proteins = protein.searchProteins(search=acc_search,
                species=selected_species, sequence=pep_search)
        
        for p in proteins:
            mods = modifications.getMeasuredPeptidesByProtein(p.id, request.user)
            peps = set()
            for mod in mods:
                [ peps.add(pep.peptide) for pep in mod.peptides ]
            protein_mods[p.id] = peps
    

    proteins = sorted(proteins, key=lambda prot: prot.acc_gene)
    
    for pid in protein_mods:
        protein_mods[pid] = sorted(protein_mods[pid], key=lambda mod: mod.getName())
        protein_mods[pid] = [ { 'site': pep.getName(), 'peptide': pep.getPeptide() } for pep in protein_mods[pid]]
    
    species_list = [ species.name for species in taxonomies.getAllSpecies() ]
    
    return {'pageTitle': strings.protein_search_page_title,
            'species_list': species_list,
            'acc_search':acc_search,
            'pep_search':pep_search,
            'stringency':stringency,
            'selected_species':selected_species,
            'proteins':proteins,
            'include_predictions': False,
            'modifications':protein_mods,
            'submitted': submitted}
