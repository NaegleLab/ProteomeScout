import celery
import logging
from ptmworker import notify_tasks
from ptmworker.helpers import upload_helpers
from ptmscout.database import modifications, protein, experiment
from ptmscout.utils import mail, uploadutils
import datetime
import traceback
log = logging.getLogger('ptmscout')

def load_scansite_peptide(pep, taxonomy):
    motif_class = None
    if 'mammalia' in taxonomy:
        motif_class="MAMMALIAN"
    elif 'saccharomycotina' in taxonomy:
        motif_class="YEAST"
    elif 'saccharomyces' in taxonomy:
        motif_class="YEAST"
    
    if motif_class != None:
        pep.predictions = upload_helpers.query_peptide_predictions(pep.pep_aligned, motif_class)
        pep.scansite_date = datetime.datetime.now()

def load_new_peptide(prot_id, pep, taxonomy):
    pep.protein_domain = protein.getProteinDomain(prot_id, pep.site_pos)
    load_scansite_peptide(pep, taxonomy)
    pep.save()


def create_errors_for_runs(exp_id, protein_accession, pep_seq, msg, runs):
    log.warning("Error importing peptide (%s, %s): %s", protein_accession, pep_seq, msg)
    for line, _rn, _s in runs:
        experiment.createExperimentError(exp_id, line, protein_accession, pep_seq, msg)


def load_protein(accession, protein_information):
    _a, _b, taxonomy, species, host_organism, _c, _d, _m, seq = protein_information
    prot = protein.getProteinBySequence(seq, species)

    if host_organism:
        taxonomy += upload_helpers.get_taxonomic_lineage(host_organism)
    
    return prot.id, prot.sequence, species, taxonomy


def load_peptide_modification(exp_id, load_ambiguities, protein_accession, protein_info, pep_seq, mods, units, series_header, runs):
    try:
        protein_id, protein_sequence, species, taxonomy = load_protein(protein_accession, protein_info)
        mod_types, aligned_sequences = upload_helpers.parse_modifications(protein_sequence, pep_seq, mods, taxonomy)

        filter_mods = []
        for i in xrange(0, len(mod_types)):
            _, seq, _ = aligned_sequences[i]
            mod = mod_types[i]
            filter_mods.append((seq, mod))

        pep_measurement = modifications.getMeasuredPeptide(exp_id, pep_seq, protein_id, filter_mods)

        if pep_measurement==None:
            pep_measurement = modifications.MeasuredPeptide()
            pep_measurement.experiment_id = exp_id
            pep_measurement.peptide = pep_seq
            pep_measurement.protein_id = protein_id

            if load_ambiguities:
                upload_helpers.check_ambiguity( pep_measurement, species )

        for i in xrange(0, len(aligned_sequences)):
            mod_type = mod_types[i] 
            site_pos, pep_sequence, _ = aligned_sequences[i]

            pep, created = upload_helpers.get_peptide(protein_id, site_pos, pep_sequence)

            if not pep_measurement.hasPeptideModification(pep, mod_type):
                pepmod = modifications.PeptideModification()
                pepmod.modification = mod_type
                pepmod.peptide = pep
                pep_measurement.peptides.append(pepmod)

            if created:
                load_new_peptide(protein_id, pep, taxonomy)

        for line, run_name, series in runs:
            upload_helpers.insert_run_data(pep_measurement, line, units, series_header, run_name, series)

        pep_measurement.save()

    except uploadutils.ParseError, pe:
        create_errors_for_runs(exp_id, protein_accession, pep_seq, pe.msg, runs)
    except TypeError:
        create_errors_for_runs(exp_id, protein_accession, pep_seq, "Failed to get data for protein: %s" % (protein_accession), runs)
    except Exception, e:
        log.warning("Unexpected Error: %s\n%s\nDuring import of peptide %d %s %s", str(e), traceback.format_exc(), exp_id, protein_accession, pep_seq)
        create_errors_for_runs(exp_id, protein_accession, pep_seq, "Unexpected error: " + str(e), runs)


UPDATE_EVERY=30

@celery.task
@upload_helpers.transaction_task
def run_peptide_import(prot_map, exp_id, peptides, mod_map, data_runs, headers, units, load_ambiguities):
    accessions = set( prot_map.keys() ) & set( peptides.keys() )

    total_peptides = 0
    for acc in accessions:
        total_peptides += len(peptides[acc])

    upload_helpers.store_stage_input(exp_id, 'peptides', prot_map)
    notify_tasks.set_loading_stage.apply_async((exp_id, 'peptides', total_peptides))

    i = 0
    for acc in accessions:
        pep_tasks = []

        for pep in peptides[acc]:
            key = (acc, pep)
            for mod_str in mod_map[key]:
                run_key = (acc, pep, mod_str)

                run_tasks = []
                for run_name in data_runs[run_key]:
                    line, series = data_runs[run_key][run_name]
                    run_tasks.append( (line, run_name, series) )

                load_peptide_modification(exp_id, load_ambiguities, acc, prot_map[acc], pep, mod_str, units, headers, run_tasks )

                i+=1
                if i % UPDATE_EVERY == 0:
                    notify_tasks.set_progress.apply_async((exp_id, i, total_peptides))

    notify_tasks.set_progress.apply_async((exp_id, i, total_peptides))
