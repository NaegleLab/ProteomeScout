from scripts.DB_init import DatabaseInitialization
from ptmscout.config import settings
from ptmscout.database import experiment, modifications, user, protein
from sqlalchemy.sql.expression import or_
from collections import defaultdict
import traceback
import os
import cPickle
import json, base64

phosphosite_exp_id = 1665 
hprd_exp_id = 1395
uniprot_exp_id = 1667 

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

def mods_by_source(mod_list):
    by_source = defaultdict(lambda: 0)
    for ms in mod_list:
        for modpep in ms.peptides:
            if ms.experiment.type=='compendia':
                by_source[ms.experiment.name] += 1
            elif ms.experiment.type=='experiment':
                by_source['experiment'] += 1
    return dd_to_dict(by_source)

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

def top_n(entries, n):
    sorted_entries = sorted(entries, key=lambda item: -item[1])

    if(len(sorted_entries) > n):
        top_n = sorted_entries[:n]
        rest_q = sum([k for (e, k) in sorted_entries[n:]])

        return top_n + [('Other', rest_q)]

    return sorted_entries

def trim_species_name(name):
    index = name.find('(')
    if index > -1:
        name = name[0:index-1]
    return name

def trim_species_names(entries):
    return [(trim_species_name(name), k) for name, k in entries]

def get_all_ptms_for_exp(ptms, eid):
    result = set()
    for ms in ptms:
        if ms.experiment_id == eid:
            for modpep in ms.peptides:
                result.add((modpep.peptide_id, modpep.modification_id))
    return result

def get_all_ptms_for_exp_type(ptms, type='compendia'):
    result = set()
    for ms in ptms:
        if ms.experiment.type == type:
            for modpep in ms.peptides:
                result.add((modpep.peptide_id, modpep.modification_id))
    return result

def compare_experiments_to_compendia(ptms):
    compendia_mods = get_all_ptms_for_exp_type(ptms, 'compendia')
    experiment_mods = get_all_ptms_for_exp_type(ptms, 'experiment')

    sets = [{'label': 'compendia', 'size': len(compendia_mods)},\
            {'label': 'experiments', 'size': len(experiment_mods)}]

    overlaps = [{'sets':[0,1], 'size':len(compendia_mods & experiment_mods)}]

    return sets, overlaps

def compare_experiments(ptms, expids, experiments):
    expids = sorted(expids)

    ptm_sets = {}
    for eid in expids:
        ptm_sets[eid] = get_all_ptms_for_exp(ptms, eid)

    sets = []
    overlaps = []

    for eid in expids:
        sets.append({'label':experiments[eid].name, 'size': len(ptm_sets[eid])})

    for i in xrange(0, len(expids)-1):
        for j in xrange(i+1, len(expids)):
            overlaps.append({'sets':[i,j], 'size':len(ptm_sets[expids[i]] & ptm_sets[expids[j]])})

    e1 = expids[0]
    e2 = expids[1]
    e3 = expids[2]
    overlaps.append({'sets':[0,1,2], 'size':len(ptm_sets[e1] & ptm_sets[e2] & ptm_sets[e3])})

    return sets, overlaps

def countMeasuredPeptides(all_mods):
    return len(all_mods)

def countReportedModifications(all_mods):
    site_set = set()
    for ms in all_mods:
        for modpep in ms.peptides:
            site_set.add((modpep.peptide_id,modpep.modification_id))
    return len(site_set)

def countSitesOfModification(all_mods):
    site_set = set()
    for ms in all_mods:
        for modpep in ms.peptides:
            site_set.add(modpep.peptide_id)
    return len(site_set)

def countProteins(all_mods):
    pidset=set()
    for ms in all_mods:
        pidset.add(ms.protein_id)
    return len(pidset)

def count_distinct_species(all_mods):
    species_ids = set()
    for ms in all_mods:
        species_ids.add(ms.protein.species_id)
    return len(species_ids)

if __name__ == "__main__":
    try:
        config_options = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        DatabaseInitialization.setUpClass(config_options)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        print "Fetching modifictions from database..."
        all_mods = dbinit.session.query(modifications.MeasuredPeptide).join(experiment.Experiment).filter(or_(experiment.Experiment.type=='compendia',experiment.Experiment.type=='experiment'),experiment.Experiment.public==1).all()

        print "Counting experiments..."
        compendia, experiments, datasets = experiment.countExperiments()
        print "Counting proteins..."
        proteins = countProteins(all_mods)
        print "Counting peptides..."
        peps = countMeasuredPeptides(all_mods)
        print "Counting modifications..."
        mods = countReportedModifications(all_mods)
        print "Counting sites..."
        sites = countSitesOfModification(all_mods)
        print "Counting users..."
        users = user.countUsers()

        print "Counting species..."
        species = count_distinct_species(all_mods)

        print "Counting modification types..."
        by_source = mods_by_source(all_mods)
        by_residue = mods_by_residue(all_mods)
        by_type = mods_by_type(all_mods)
        by_species = mods_by_species(all_mods)
        ptms = count_distinct_ptm_types(all_mods)

        print "Formatting PTMs"
        all_ptms = modifications.getAllPTMs()
        formatted_ptms = format_ptms(all_ptms)

        by_residue_short = top_n(by_residue.items(), 5)
        by_source_short = by_source.items()
        by_type_short = top_n(by_type.items(), 5)
        by_species_short = trim_species_names( top_n(by_species.items(), 5) )


        experiments_dict = {}
        for exp in dbinit.session.query(experiment.Experiment).filter(experiment.Experiment.public == 1):
            experiments_dict[exp.id] = exp

        sets, overlaps = compare_experiments( all_mods, [ hprd_exp_id, phosphosite_exp_id, uniprot_exp_id ], experiments_dict )
        compendia_venn_diagram = {'sets':sets, 'overlaps': overlaps}

        sets, overlaps = compare_experiments_to_compendia(all_mods)
        experiment_venn_diagram = {'sets':sets, 'overlaps': overlaps}

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
                'ptm_defs': formatted_ptms,
                'by_residue_json': base64.b64encode(json.dumps(by_residue_short)),
                'by_source_json': base64.b64encode(json.dumps(by_source_short)),
                'by_type_json': base64.b64encode(json.dumps(by_type_short)),
                'by_species_json': base64.b64encode(json.dumps(by_species_short)),
                'mod_venn_diagram': base64.b64encode(json.dumps(compendia_venn_diagram)),
                'exp_venn_diagram': base64.b64encode(json.dumps(experiment_venn_diagram)),
                }

        with open( os.path.join(settings.ptmscout_path, settings.statistics_file), 'w' ) as pypfile:
            cPickle.dump(rval, pypfile)

    except Exception, e:
        traceback.print_exc()
        dbinit.rollback()
    else:
        dbinit.tearDown()
