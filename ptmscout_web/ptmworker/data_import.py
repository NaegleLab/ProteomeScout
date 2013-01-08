import celery
import logging
from ptmworker.helpers import upload_helpers
from ptmworker import notify_tasks, protein_tasks, GO_tasks, peptide_tasks
from ptmscout.database import upload, modifications, protein, experiment

log = logging.getLogger('ptmscout')

@celery.task
@upload_helpers.transaction_task
def start_import(exp_id, session_id, user_email, application_url):
    exp = experiment.getExperimentById(exp_id, check_ready=False, secure=False)

    log.info("Loading session info...")
    session = upload.getSessionById(session_id, secure=False)
    
    log.info("Loading data file...")
    accessions, peptides, mod_map, data_runs, errors, line_mapping = upload_helpers.parse_datafile(session)

    if exp.loading_stage == 'query' or exp.status != 'error':
        exp.clearErrors()
        log.info("Reporting data file errors...")
        upload_helpers.report_errors(exp_id, errors, line_mapping)
        exp.saveExperiment()

    if len(accessions) == 0:
        log.info("Nothing to do: all proteins were rejected!")
        notify_tasks.finalize_import.apply_async((exp_id, user_email, application_url))
    else:
        headers = upload_helpers.get_series_headers(session)

        notify_tasks.set_loading_status.apply_async((exp_id, 'loading'))

        query_task = protein_tasks.get_proteins_from_external_databases.s(accessions, exp_id, line_mapping)
        proteins_task = protein_tasks.query_protein_metadata.s(accessions, exp_id, line_mapping)
        GO_task = GO_tasks.import_go_terms.s(exp_id)
        peptide_task = peptide_tasks.run_peptide_import.s(exp_id, peptides, mod_map, data_runs, headers, session.units)
        finalize_task = notify_tasks.finalize_import.si(exp_id, user_email, application_url)

        last_stage_arg = exp.get_last_result()

        if exp.loading_stage == 'query':
            load_task = ( query_task | proteins_task | GO_task | peptide_task | finalize_task )
        elif exp.loading_stage == 'proteins':
            load_task = ( proteins_task | GO_task | peptide_task | finalize_task )
        elif exp.loading_stage == 'GO terms':
            load_task = ( GO_task | peptide_task | finalize_task )
        elif exp.loading_stage == 'peptides':
            load_task = ( peptide_task | finalize_task )

        load_task.apply_async((last_stage_arg,), link_error=notify_tasks.finalize_experiment_error_state.s(exp_id, user_email, application_url))

        log.info("Tasks started... now we wait")
    
