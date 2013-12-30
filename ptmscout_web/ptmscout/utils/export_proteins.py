from collections import defaultdict
from ptmscout.utils import protein_utils
from ptmscout.config import settings, strings

def format_protein_accessions(accessions, query_accessions):
    accessions = [ acc.value for acc in sorted(accessions, key=lambda acc: acc.value) if protein_utils.get_accession_type( acc.value ) in protein_utils.get_valid_accession_types() ]
    accessions = sorted(accessions, key=lambda item: (0 if item in query_accessions else 1, item))
    return "; ".join(accessions)

def check_modtype_filter(mod, modtype_filter):
    if modtype_filter == None:
        return True

    modtype_filter = modtype_filter.lower()

    p = mod
    while(p is not None):
        if p.name.lower() == modtype_filter:
            return True
        for k in p.keywords:
            if k.keyword.lower() == modtype_filter:
                return True
        p = p.parent

    return False

def get_query_accessions(mods):
    return set([ ms.query_accession for ms in mods ])

def filter_sites(ms, regions):
    return list( set([ region
                    for modpep in ms.peptides
                    for region in regions
                    if region.hasSite(modpep.peptide.site_pos) ]))

def filter_site_regions(ms, regions, types):
    return list( set([ region
                    for modpep in ms.peptides
                    for region in regions
                    if region.hasSite(modpep.peptide.site_pos) and region.type in types ]))

def filter_regions(regions, types):
    return list( set([ region
                    for region in regions
                    if region.type in types ]))


def format_modifications(mods, modtype_filter):
    modlist = [ (ms.experiment_id, modpep.peptide, modpep.modification) for ms in mods for modpep in ms.peptides ]
    modlist = [ (exp_id, pep.site_pos, "%s-%s" % (pep.getName(), mod.name)) for exp_id, pep, mod in modlist if check_modtype_filter(mod, modtype_filter) ]

    explist = defaultdict(set)
    for exp_id, site_pos, modstr in modlist:
        explist[(site_pos, modstr)].add(exp_id)

    modlist = [ (site_pos, modstr) for _, site_pos, modstr in modlist ]
    modlist = [ modstr for site_pos, modstr in sorted( list( set(modlist) ) ) ]

    explist = [ ','.join([ str(exp_id) for exp_id in sorted( list( explist[k] )) ]) for k in sorted(explist.keys())  ]

    n = len(modlist)
    return n, '; '.join(modlist), '; '.join(explist)

def format_regions(regions):
    return (settings.mod_separator_character + ' ').join( [ "%s:%d-%d" % (r.label, r.start, r.stop) for r in sorted( regions, key=lambda r: r.start ) ] )

def format_region_types(regions):
    return (settings.mod_separator_character + ' ').join( [ "%s:%d-%d" % (r.type, r.start, r.stop) for r in sorted( regions, key=lambda r: r.start ) ] )

def format_mutations(mutations):
    return (settings.mod_separator_character + ' ').join( [ str(m) for m in sorted(mutations, key=lambda m: m.location) ] )

def format_GO_terms(prot):
    return (settings.mod_separator_character + ' ').join( [ "%s-%s" % (goe.GO_term.GO, goe.GO_term.term) for goe in sorted(prot.GO_terms, key=lambda term: term.GO_term.GO) ] )

def format_scansite(mods):
    plist = [ (modpep.peptide.site_pos, modpep.peptide.site_type, pred.source, pred.value, pred.percentile) for ms in mods for modpep in ms.peptides for pred in modpep.peptide.predictions ]
    return '; '.join( [ '%s%d-%s-%s:%.2f' % (site_tp, pos, source, value, percentile) for pos, site_tp, source, value, percentile in sorted( list( set( plist ) ) ) ])


