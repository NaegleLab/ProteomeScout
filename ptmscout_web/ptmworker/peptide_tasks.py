import celery
import logging
from ptmworker import notify_tasks
from ptmworker.helpers import upload_helpers, scansite_tools
from ptmscout.database import modifications, protein, experiment
from ptmscout.utils import uploadutils
import datetime
import traceback
from sqlalchemy.exc import DBAPIError, SQLAlchemyError

log = logging.getLogger('ptmscout')

def load_scansite_peptide(pep, taxonomy):
    scansite_residues = set(['Y', 'S', 'T'])
    motif_class = None
    if 'mammalia' in taxonomy:
        motif_class="MAMMALIAN"
    elif 'saccharomycotina' in taxonomy:
        motif_class="YEAST"
    elif 'saccharomyces' in taxonomy:
        motif_class="YEAST"
    
    if motif_class != None and pep.site_type in scansite_residues:
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

def load_protein(accession, protein_record):
    prot = protein.getProteinBySequence(protein_record.sequence, protein_record.species)
    return prot.id, prot.sequence, protein_record.species, protein_record.full_taxonomy



def load_peptide_modification(exp_id, load_ambiguities, protein_accession, protein_record, site_designation, mods, nullmods, units, series_header, runs, is_site=False):
    try:
        protein_id, protein_sequence, species, taxonomy = load_protein(protein_accession, protein_record)

        pep_seq = site_designation
        if is_site:
            pep_seq = upload_helpers.get_pep_seq_from_sites(protein_sequence, site_designation)

        mod_types = []
        aligned_sequences = []
        if nullmods:
            aligned_sequences = upload_helpers.parse_nullmod(protein_sequence, pep_seq)
            none_mod = modifications.getModificationByName('None')
            mod_types = [ none_mod for _seq in aligned_sequences ]
        else:
            mod_types, aligned_sequences = upload_helpers.parse_modifications(protein_sequence, pep_seq, mods, taxonomy)

        filter_mods = []
        for i in xrange(0, len(mod_types)):
            _, seq, _ = aligned_sequences[i]
            mod = mod_types[i]
            filter_mods.append((seq, mod))

        pep_measurement = modifications.getMeasuredPeptide(exp_id, site_designation, protein_id, filter_mods)

        if pep_measurement==None:
            pep_measurement = modifications.MeasuredPeptide()
            pep_measurement.experiment_id = exp_id
            pep_measurement.peptide = site_designation
            pep_measurement.protein_id = protein_id
            pep_measurement.query_accession = protein_accession

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
        create_errors_for_runs(exp_id, protein_accession, site_designation, pe.msg, runs)
    except scansite_tools.ScansiteError:
        create_errors_for_runs(exp_id, protein_accession, site_designation, "Scansite query failed for peptide: %s '%s'" % (protein_accession, site_designation), runs)
    except TypeError:
        create_errors_for_runs(exp_id, protein_accession, site_designation, "Failed to get data for protein: %s" % (protein_accession), runs)
    except DBAPIError:
        raise
    except SQLAlchemyError:
        raise
    except Exception, e:
        log.warning("Unexpected Error: %s\n%s\nDuring import of peptide %d %s %s", str(e), traceback.format_exc(), exp_id, protein_accession, site_designation)
        create_errors_for_runs(exp_id, protein_accession, site_designation, "Unexpected error: " + str(e), runs)


UPDATE_EVERY=30

@celery.task
@upload_helpers.notify_job_failed
@upload_helpers.transaction_task
def run_peptide_import(prot_map, peptides, mod_map, data_runs, headers, units, load_ambiguities, nullmods, by_site, exp_id, job_id):
    accessions = set( prot_map.keys() ) & set( peptides.keys() )

    total_peptides = 0
    for acc in accessions:
        total_peptides += len(peptides[acc])

    upload_helpers.store_stage_input(exp_id, 'peptides', prot_map)
    notify_tasks.set_job_stage.apply_async((job_id, 'peptides', total_peptides))

    i = 0
    for acc in accessions:
        for pep in peptides[acc]:
            key = (acc, pep)
            for mod_str in mod_map[key]:
                run_key = (acc, pep, mod_str)

                run_tasks = []
                for run_name in data_runs[run_key]:
                    line, series = data_runs[run_key][run_name]
                    run_tasks.append( (line, run_name, series) )

                load_peptide_modification(exp_id, load_ambiguities and not by_site, acc, prot_map[acc], pep, mod_str, nullmods, units, headers, run_tasks, is_site = by_site)

                i+=1
                if i % UPDATE_EVERY == 0:
                    notify_tasks.set_job_progress.apply_async((job_id, i, total_peptides))

    notify_tasks.set_job_progress.apply_async((job_id, i, total_peptides))
