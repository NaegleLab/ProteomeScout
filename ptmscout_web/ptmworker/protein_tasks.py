import celery
import logging 
from ptmworker import notify_tasks
from ptmworker.helpers import upload_helpers, entrez_tools, pfam_tools, quickgo_tools, picr_tools, uniprot_tools
from ptmscout.database import upload, modifications, protein, experiment
from celery.canvas import group
from ptmscout.config import strings, settings
from ptmscout.utils import mail, uploadutils
import datetime
import traceback
from ptmscout.utils.decorators import rate_limit

log = logging.getLogger('ptmscout')

MAX_NCBI_BATCH_SIZE = 500
MAX_UNIPROT_BATCH_SIZE = 200

def get_uniprot_proteins(protein_accessions):
    log.info("Getting uniprot records for %d accessions", len(protein_accessions))
    prot_map = uniprot_tools.get_uniprot_records(protein_accessions)

    errors = []
    for acc in protein_accessions:
        if acc not in prot_map:
            errors.append(entrez_tools.EntrezError())
            errors[-1].acc = acc

    return prot_map, errors

def log_errors(query_errors, exp_id, accessions, line_mappings):
    log.info("Detected %d errors", len(query_errors))
    
    #report errors
    for error in query_errors:
        for line in accessions[error.acc]:
            accession, peptide = line_mappings[line]
            experiment.createExperimentError(exp_id, line, accession, peptide, strings.experiment_upload_warning_accession_not_found % (accession))


def load_new_protein(accession, protein_information):
    name, gene, taxonomy, species, host_organism, prot_accessions, domains, seq = protein_information

    prot = protein.getProteinBySequence(seq, species)
    if prot == None:
        prot = upload_helpers.create_new_protein(name, gene, seq, species)

    # load additional protein accessions if available

    other_accessions = picr_tools.get_picr(accession, prot.species.taxon_id)
    added_accessions = upload_helpers.create_accession_for_protein(prot, prot_accessions + other_accessions)

    # viruses need the taxonomy of their host organism to query scansite and
    # check for valid PTMs

    if host_organism:
        taxonomy += upload_helpers.get_taxonomic_lineage(host_organism)

    pfam_tools.parse_or_query_domains(prot, domains, accession)
    upload_helpers.map_expression_probesets(prot)

    prot.saveProtein()

    log.info("Created protein: %s | %s" , accession, str(added_accessions))
    return prot

def report_protein_error(acc, protein_map, accessions, line_mappings, exp_id, message):
    log.warning("Failed to create protein: %s -- '%s'", acc, message)
    del protein_map[acc]

    if acc in accessions:
        for line in accessions[acc]:
            accession, peptide = line_mappings[line]
            experiment.createExperimentError(exp_id, line, accession, peptide, message)

UPDATE_EVERY = 100
def create_missing_proteins(protein_map, missing_proteins, accessions, exp_id, line_mappings):

    i = 0
    #create entries for the missing proteins
    protein_id_map = {}
    for acc in missing_proteins:
        try:
            prot = load_new_protein(acc, protein_map[acc])
            protein_id_map[acc] = prot.id

        except uploadutils.ParseError, e:
            report_protein_error(acc, protein_map, accessions, line_mappings, exp_id, e.message)
        except pfam_tools.PFamError:
            report_protein_error(acc, protein_map, accessions, line_mappings, exp_id, "PFam query failed for protein: %s" % (acc))
        except picr_tools.PICRError:
            report_protein_error(acc, protein_map, accessions, line_mappings, exp_id, "PICR query failed for protein: %s" % (acc))
        except Exception, e:
            log.warning("Unexpected Error: %s\n%s\nDuring import of protein %s", str(e), traceback.format_exc(), acc)
            report_protein_error(acc, protein_map, accessions, line_mappings, exp_id, "Unexpected Error: %s" % (str(e)))

        i+=1
        if i % UPDATE_EVERY == 0:
            notify_tasks.set_progress.apply_async((exp_id, i, len(missing_proteins)))

    notify_tasks.set_progress.apply_async((exp_id, i, len(missing_proteins)))
    return protein_map, protein_id_map


@celery.task
@upload_helpers.transaction_task
def query_protein_metadata(external_db_result, accessions, exp_id, line_mapping):
    #list the missing proteins
    missing_proteins = set()
    for acc in external_db_result:
        _n, _g, _t, species, _h, _a, _d, seq = external_db_result[acc]
        if protein.getProteinBySequence(seq, species) == None:
            missing_proteins.add(acc)

    notify_tasks.set_loading_stage.apply_async((exp_id, 'proteins', external_db_result, len(missing_proteins)))
    return create_missing_proteins(external_db_result, missing_proteins, accessions, exp_id, line_mapping)

@celery.task
@upload_helpers.transaction_task
def get_proteins_from_external_databases(ignored, accessions, exp_id, line_mapping):
    uniprot_ids, other_ids = upload_helpers.extract_uniprot_accessions(accessions.keys())

    uniprot_tasks = upload_helpers.create_chunked_tasks_preserve_groups(sorted(uniprot_ids), MAX_UNIPROT_BATCH_SIZE)
    ncbi_tasks = upload_helpers.create_chunked_tasks(sorted(other_ids), MAX_NCBI_BATCH_SIZE)
    total_task_cnt = len(uniprot_tasks) + len(ncbi_tasks)

    notify_tasks.set_loading_stage.apply_async((exp_id, 'query', None, total_task_cnt))

    i = 0
    protein_map = {}
    for ncbi_accessions in ncbi_tasks:
        log.info("Getting Geeneus records for %d accessions", len(ncbi_accessions))
        result, errors = entrez_tools.get_proteins_from_ncbi(ncbi_accessions)
        protein_map.update(result)

        log_errors(errors, exp_id, accessions, line_mapping)

        i+=1
        notify_tasks.set_progress.apply_async((exp_id, i, total_task_cnt))

    for uniprot_accessions in uniprot_tasks:
        result, errors = get_uniprot_proteins(uniprot_accessions)
        protein_map.update(result)

        log_errors(errors, exp_id, accessions, line_mapping)

        i+=1
        notify_tasks.set_progress.apply_async((exp_id, i, total_task_cnt))

    return protein_map
