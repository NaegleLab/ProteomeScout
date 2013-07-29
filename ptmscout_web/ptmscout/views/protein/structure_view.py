from pyramid.view import view_config
from ptmscout.config import strings, settings
from ptmscout.database import protein, modifications
from ptmscout.views.protein import decorators
import json, base64
from collections import defaultdict

def format_scansite_predictions(prot):
    formatted_scansite = defaultdict(list)
    
    for pred in prot.scansite:
        ss = { 'source': pred.source,
               'value': pred.value,
               'score': "%.2f%%" % (pred.percentile)
               }
        formatted_scansite[pred.site_pos].append( ss )

    return formatted_scansite

def format_protein_mutations(prot):
    formatted_mutations = {}
    for m in prot.mutations:
        mut_dict = {}
        mut_dict['type'] = m.mutationType
        mut_dict['location'] = m.location
        mut_dict['original'] = m.original
        mut_dict['mutant'] = m.mutant
        mut_dict['annotation'] = m.annotation

        mut_list = formatted_mutations.get(m.location, [])
        mut_list.append(mut_dict)
        formatted_mutations[m.location] = mut_list

    return formatted_mutations

def get_activation_loops(prot):
    formatted_regions = []
    for d in prot.regions:
        if d.type == 'Activation Loop':
            region_dict = {}
            region_dict['type'] = d.type
            region_dict['label'] = d.label
            region_dict['source'] = d.source
            region_dict['start'] = d.start
            region_dict['stop'] = d.stop

            formatted_regions.append(region_dict)

    return sorted(formatted_regions, key=lambda d: d['start'])

def get_uniprot_domains(prot):
    formatted_regions = []

    for d in prot.regions:
        if d.type == 'domain' and d.source == 'uniprot':
            region_dict = {}
            region_dict['type'] = d.type
            region_dict['label'] = d.label
            region_dict['source'] = d.source
            region_dict['start'] = d.start
            region_dict['stop'] = d.stop

            formatted_regions.append(region_dict)

    return sorted(formatted_regions, key=lambda d: d['start'])

def get_ncbi_domains(prot):
    formatted_regions = []

    for d in prot.regions:
        if d.type == 'Domain' and d.source == 'ncbi':
            region_dict = {}
            region_dict['type'] = d.type
            region_dict['label'] = d.label
            region_dict['source'] = d.source
            region_dict['start'] = d.start
            region_dict['stop'] = d.stop

            formatted_regions.append(region_dict)

    return sorted(formatted_regions, key=lambda d: d['start'])

def format_protein_regions(prot):
    formatted_regions = {}

    formatted_regions['activation_loops'] = get_activation_loops(prot)
    formatted_regions['uniprot_domains'] = get_uniprot_domains(prot)
    formatted_regions['ncbi_domains'] = get_ncbi_domains(prot)

    return formatted_regions

def format_protein_domains(prot):
    formatted_domains = []
    for d in prot.domains:
        domain_dict = {}
        domain_dict['label'] = d.label
        domain_dict['source'] = d.source
        domain_dict['start'] = d.start
        domain_dict['stop'] = d.stop

        formatted_domains.append(domain_dict)

    return sorted(formatted_domains, key=lambda d: d['start'])

def get_site_regions(regions, pos):
    region_names = []
    for r in regions:
        if r.start <= pos and pos <= r.stop:
            region_names.append(r.label)

    return region_names

def format_protein_modifications(request, prot, mod_sites):
    experiments = {}
    mod_types = set()
    mods = {}

    for ms in mod_sites:
        experiments[ms.experiment_id] = ms.experiment.name

        exp_url = request.route_url('experiment', id=ms.experiment_id)
        if ms.experiment.type == 'compendia':
            exp_url = ms.experiment.URL

        for modpep in ms.peptides:
            mod_types.add(modpep.modification.name)
            pep = modpep.peptide
            ptm = modpep.modification

            pos_description = mods.get( pep.site_pos,
                                        { 'mods': {},
                                          'residue': pep.site_type,
                                          'domain': pep.protein_domain.label if pep.protein_domain else None,
                                          'regions': get_site_regions(prot.regions, pep.site_pos),
                                          'peptide': pep.pep_aligned} )

            mod_list = pos_description['mods'].get(ptm.name, [])
            mod_list.append( {'MS': ms.id, 'experiment_url': exp_url, 'experiment': ms.experiment_id, 'has_data': len(ms.data) > 0 } )
            pos_description['mods'][ptm.name] = mod_list

            mods[pep.site_pos] = pos_description

    return experiments, sorted(list(mod_types)), mods


@view_config(route_name='protein_viewer', renderer='ptmscout:templates/proteins/protein_viewer.pt')
@decorators.experiment_filter
def protein_structure_viewer(context, request):
    protein_id = int(request.matchdict['id'])
    prot = protein.getProteinById(protein_id)

    mod_sites = modifications.getMeasuredPeptidesByProtein(protein_id, request.user)

    formatted_exps, formatted_mod_types, formatted_mods = format_protein_modifications(request, prot, mod_sites)
    formatted_domains = format_protein_domains(prot)
    formatted_regions = format_protein_regions(prot)
    formatted_mutations = format_protein_mutations(prot)
    formatted_scansite = format_scansite_predictions(prot)

    data = {'seq': prot.sequence, 
            'domains': formatted_domains, 
            'mods': formatted_mods,
            'mutations': formatted_mutations,
            'regions': formatted_regions,
            'mod_types': formatted_mod_types,
            'scansite': formatted_scansite,
            'exps': formatted_exps,
            'pfam_url': settings.pfam_family_url,
            'protein_data_url': "%s/proteins/%d/data" % (request.application_url, protein_id),
            'experiment': request.urlfilter.get_field('experiment_id') }
    encoded_data = base64.b64encode( json.dumps( data ) )

    return {'pageTitle': strings.protein_structure_page_title,
            'protein': prot,
            'experiments': formatted_exps,
            'mod_types': formatted_mod_types,
            'tracks': ["PFam Domains", "PTMs", "Activation Loops", "Uniprot Domains", "Entrez Domains", "Mutations", "Scansite"],
            'data':encoded_data}
