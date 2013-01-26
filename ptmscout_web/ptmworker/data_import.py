import celery
import logging
from ptmworker.helpers import upload_helpers
from ptmworker import notify_tasks, protein_tasks, GO_tasks, peptide_tasks
from ptmscout.database import upload, modifications, protein, experiment
import transaction
log = logging.getLogger('ptmscout')

@celery.task
@upload_helpers.dynamic_transaction_task
def start_import(exp_id, session_id, user_email, application_url):
    exp = experiment.getExperimentById(exp_id, check_ready=False, secure=False)

    log.info("Loading session info...")
    session = upload.getSessionById(session_id, secure=False)
    
    log.info("Loading data file...")
    accessions, peptides, mod_map, data_runs, errors, line_mapping = upload_helpers.parse_datafile(session)

    if exp.loading_stage == 'in queue' or exp.status != 'error':
        #exp.clearErrors()
        log.info("Reporting data file errors...")
        upload_helpers.report_errors(exp_id, errors, line_mapping)
        exp.saveExperiment()

    if len(accessions) == 0:
        log.info("Nothing to do: all proteins were rejected!")
        notify_tasks.finalize_import.apply_async((exp_id, user_email, application_url))
    else:
        headers = upload_helpers.get_series_headers(session)

        exp.status = 'loading'
        experiment.setExperimentProgress(exp_id, 0, 0)
        exp.saveExperiment()

        load_ambiguities = exp.ambiguity == 1

        query_task = protein_tasks.get_proteins_from_external_databases.s(accessions, exp_id, line_mapping)
        proteins_task = protein_tasks.query_protein_metadata.s(accessions, exp_id, line_mapping)
        GO_task = GO_tasks.import_go_terms.s(exp_id)
        peptide_task = peptide_tasks.run_peptide_import.s(exp_id, peptides, mod_map, data_runs, headers, session.units, load_ambiguities)
        finalize_task = notify_tasks.finalize_import.si(exp_id, user_email, application_url)

        last_stage_arg = upload_helpers.get_stage_input(exp.id, exp.loading_stage)

        if exp.loading_stage == 'query':
            load_task = ( query_task | proteins_task | GO_task | peptide_task | finalize_task )
        elif exp.loading_stage == 'proteins':
            load_task = ( proteins_task | GO_task | peptide_task | finalize_task )
        elif exp.loading_stage == 'GO terms':
            load_task = ( GO_task | peptide_task | finalize_task )
        elif exp.loading_stage == 'peptides':
            load_task = ( peptide_task | finalize_task )

        callback_task = notify_tasks.finalize_experiment_error_state.s(exp_id, user_email, application_url)

        log.info("Tasks created... now we wait")

        return load_task, last_stage_arg, callback_task
