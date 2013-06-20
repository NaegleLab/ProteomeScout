from scripts.DB_init import DatabaseInitialization
from ptmscout.config import settings
from ptmscout.database import experiment, modifications, user, protein
from collections import defaultdict
import traceback
import os
import cPickle

def get_sorted_children(parent_id, ptms):
    children = [ p for p in ptms if p.parent_id == parent_id ]
    ordered_ptms = []

    for c in children:
        ordered_ptms.append(c)
        ordered_ptms += get_sorted_children(c.id, ptms)

    return ordered_ptms

def format_ptms(ptms):
    ordered_ptms = get_sorted_children(None, ptms)
    formatted_ptms = []

    for ptm in ordered_ptms:
        formatted = {}
        formatted['level'] = ptm.getLevel()
        formatted['name'] = ptm.name
        formatted['position'] = ptm.position
        target_set = ptm.getTargets()
        if None in target_set:
            target_set.remove(None)
        formatted['residues'] = ', '.join( sorted(target_set) )
        formatted['species'] = ', '.join( sorted(ptm.getTaxons()) )
        formatted['aliases'] = ', '.join( sorted( [k.keyword for k in ptm.keywords] ) )
        formatted['monomass'] = ptm.mono_mass_diff
        formatted['avgmass'] = ptm.avg_mass_diff
        formatted_ptms.append(formatted)

    return formatted_ptms

def dd_to_dict(dd):
    d = dict()
    for k in dd:
        d[k] = dd[k]
    return d

def mods_by_residue(mod_list):
    by_residue = defaultdict(lambda: 0)
    for ms in mod_list:
        for modpep in ms.peptides:
            by_residue[modpep.peptide.site_type] += 1
    return dd_to_dict(by_residue)

def mods_by_type(mod_list):
    by_type = defaultdict(lambda: 0)
    for ms in mod_list:
        for modpep in ms.peptides:
            by_type[modpep.modification.name] += 1
    return dd_to_dict(by_type)

def mods_by_species(mod_list):
    by_species = defaultdict(lambda: 0)
    for ms in mod_list:
        by_species[ms.protein.species.name] += len(ms.peptides)
    return dd_to_dict(by_species)

def count_distinct_ptm_types(mod_list):
    mod_ids = set()
    for ms in mod_list:
        for modpep in ms.peptides:
            mod_ids.add(modpep.modification_id)

    return len(mod_ids)

def count_distinct_species(protein_list):
    species_ids = set()
    for pr in protein_list:
        species_ids.add(pr.species_id)
    return len(species_ids)


if __name__ == "__main__":
    try:
        config_options = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        DatabaseInitialization.setUpClass(config_options)
        dbinit = DatabaseInitialization()
        dbinit.setUp()


        print "Counting experiments..."
        compendia, experiments, datasets = experiment.countExperiments()
        print "Counting proteins..."
        proteins = protein.countProteins()
        print "Counting peptides..."
        peps = modifications.countMeasuredPeptides()
        print "Counting modifications..."
        mods = modifications.countReportedModifications()
        print "Counting sites..."
        sites = modifications.countSitesOfModification()
        print "Counting users..."
        users = user.countUsers()

        print "Counting species..."
        all_proteins = protein.getAllProteins()
        species = count_distinct_species(all_proteins)

        print "Counting modification types..."
        all_mods = modifications.getAllMeasuredPeptides()
        by_residue = mods_by_residue(all_mods)
        by_type = mods_by_type(all_mods)
        by_species = mods_by_species(all_mods)
        ptms = count_distinct_ptm_types(all_mods)

        print "Formatting PTMs"
        all_ptms = modifications.getAllPTMs()
        formatted_ptms = format_ptms(all_ptms)

        rval = {
                'compendia': compendia,
                'experiments': experiments,
                'datasets': datasets,
                'users': users,
                'proteins': proteins,
                'peps': peps,
                'mods': mods,
                'sites': sites,
                'by_residue': by_residue,
                'by_species': by_species,
                'by_type': by_type,
                'ptms': ptms,
                'species': species,
                'ptm_defs': formatted_ptms
                }

        with open( os.path.join(settings.ptmscout_path, settings.statistics_file), 'w' ) as pypfile:
            cPickle.dump(rval, pypfile)

    except Exception, e:
        traceback.print_exc()
        dbinit.rollback()
    else:
        dbinit.tearDown()
