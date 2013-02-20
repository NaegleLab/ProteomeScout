from pyramid.view import view_config
from ptmscout.config import strings, settings
from ptmscout.database import protein, modifications
from ptmscout.views.protein import decorators
import json, base64


def format_scansite_predictions(prot, mod_sites):
    peptides = {}
    for ms in mod_sites:
        for modpep in ms.peptides:
            pep = modpep.peptide
            peptides[pep.site_pos] = pep

    formatted_scansite = {}
    for site in peptides:
        pep = peptides[site]
        scansite_predictions = []
        for ss in pep.predictions:
            scansite_predictions.append({ \
                    'source': ss.source,
                    'value': ss.value,
                    'score': ss.score
                })
        formatted_scansite[site] = scansite_predictions

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


def format_protein_regions(prot):
    formatted_regions = []
    for d in prot.regions:
        region_dict = {}
        region_dict['label'] = d.label
        region_dict['source'] = d.source
        region_dict['start'] = d.start
        region_dict['stop'] = d.stop

        formatted_regions.append(region_dict)

    return sorted(formatted_regions, key=lambda d: d['start'])

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

def format_protein_modifications(prot, mod_sites):
    experiments = {}
    mod_types = set()
    mods = {}

    for ms in mod_sites:
        experiments[ms.experiment_id] = ms.experiment.name

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
            mod_list.append( {'MS': ms.id, 'experiment': ms.experiment_id, 'exported': ms.experiment.export == 1, 'has_data': len(ms.data) > 0 } )
            pos_description['mods'][ptm.name] = mod_list

            mods[pep.site_pos] = pos_description

    return experiments, sorted(list(mod_types)), mods


@view_config(route_name='protein_viewer', renderer='ptmscout:templates/proteins/protein_viewer.pt')
@decorators.experiment_filter
def protein_structure_viewer(request):
    protein_id = int(request.matchdict['id'])
    prot = protein.getProteinById(protein_id)

    mod_sites = modifications.getMeasuredPeptidesByProtein(protein_id, request.user)

    formatted_exps, formatted_mod_types, formatted_mods = format_protein_modifications(prot, mod_sites)
    formatted_domains = format_protein_domains(prot)
    formatted_regions = format_protein_regions(prot)
    formatted_mutations = format_protein_mutations(prot)
    formatted_scansite = format_scansite_predictions(prot, mod_sites)

    data = {'seq': prot.sequence, 
            'domains': formatted_domains, 
            'mods': formatted_mods,
            'mutations': formatted_mutations,
            'regions': formatted_regions,
            'mod_types': formatted_mod_types,
            'scansite': formatted_scansite,
            'exps': formatted_exps,
            'pfam_url': settings.pfam_family_url,
            'experiment_url': "%s/experiments" % (request.application_url),
            'protein_data_url': "%s/proteins/%d/data" % (request.application_url, protein_id),
            'experiment': request.urlfilter.get_field('experiment_id') }
    encoded_data = base64.b64encode( json.dumps( data ) )

    return {'pageTitle': strings.protein_structure_page_title,
            'protein': prot,
            'experiments': formatted_exps,
            'mod_types': formatted_mod_types,
            'tracks': ["Domains", "PTMs", "Regions", "Mutations"],
            'data':encoded_data}
