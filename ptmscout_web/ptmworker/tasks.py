import celery
from ptmscout.config import strings
from celery.canvas import group
from ptmscout.utils import mail, uploadutils
from ptmworker import upload_helpers
import logging 
import traceback
from ptmscout.database import modifications, experiment, upload
import transaction

log = logging.getLogger('ptmscout')

@celery.task
def finalize_import(exp_id, user_email, application_url):
    try:
        log.debug("finalization...")
        
        exp = upload_helpers.mark_experiment(exp_id, 'loaded')
        
        measuredPeps = modifications.getMeasuredPeptidesByExperiment(exp.id, check_ready=False, secure=False)
        prot_ids = set([m.protein_id for m in measuredPeps])

        for error in exp.errors:
            log.debug("Line %d error: %s", error.line, error.message)
        
        message = strings.experiment_upload_finished_message % (exp.name, len(measuredPeps), len(prot_ids), len(exp.errors), application_url + "/experiments/%d/errors" % (exp_id))
        mail.celery_send_mail([user_email], strings.experiment_upload_finished_subject, message)
    except:
        log.debug(traceback.format_exc())
        log.debug("Error in finalization")
        upload_helpers.mark_experiment(exp_id, 'error')
    
    

    
def logErrorsToDB(exp_id, affected_lines, e):
    for line in affected_lines:
        experiment.createExperimentError(exp_id, line, e.msg)
    
    
    
@celery.task
def load_peptide((protein_id, prot_seq, taxonomy), exp_id, pep_seq, modlist, run_task_args):
    affected_lines = [line for line, _, _, _, _ in run_task_args]
    
    try:
        log.debug("loading peptide: %s %s", pep_seq, modlist)
        created_mods = upload_helpers.create_modifications(protein_id, prot_seq, pep_seq, modlist, taxonomy)
        
        MS_pep = modifications.MeasuredPeptide()
        MS_pep.experiment_id = exp_id
        MS_pep.phosphopep = pep_seq
        MS_pep.protein_id = protein_id
        
        for modified_peptide in created_mods:
            MS_pep.phosphopeps.append(modified_peptide)
            
        MS_pep.save()
        
        for line, units, headers, name, series in run_task_args:
            upload_helpers.insert_run_data(MS_pep, line, units, headers, name, series)
            
    except uploadutils.ParseError, e:
        logErrorsToDB(exp_id, affected_lines, e)
    except Exception, e:
        log.debug(traceback.format_exc())
    

@celery.task
def load_protein(protein_information, exp_id, affected_lines):
    try:
        name, gene, taxonomy, species, prot_accessions, seq = protein_information
        prot = upload_helpers.find_protein(name, gene, seq, prot_accessions, species)
        
        # execute extra protein data import tasks (Future)

        protein_id = prot.id
        protein_seq = prot.sequence
        
            
        return protein_id, protein_seq, taxonomy
    except uploadutils.ParseError, e:
        logErrorsToDB(exp_id, affected_lines, e)
    except Exception, e:
        log.debug(traceback.format_exc())
    
    
@celery.task
def process_error_state(exp_id):
    exp = upload_helpers.mark_experiment(exp_id, 'error')
    log.debug("an error occurred during processing '%s'", exp.name)
    

def create_import_tasks(exp_id, prot_map, accessions, peptides, mod_map, series_headers, units, data_runs):
    import_tasks = []

    for acc in prot_map:
        pep_tasks = []
        
        for pep in peptides[acc]:
            run_tasks = []
            
            runs = data_runs[(acc, pep)]
            for run_name in runs:
                line, series =  runs[run_name]
                run_tasks.append((line, units, series_headers, run_name, series))
            
            pep_tasks.append( load_peptide.s(exp_id, pep, mod_map[(acc,pep)], run_tasks) )
        
        import_tasks.append(( load_protein.s(prot_map[acc], exp_id, accessions[acc]) | group(pep_tasks) ))

    return import_tasks

def invoke(import_tasks, exp_id, user_email, application_url):
    return ( group(import_tasks) | finalize_import.si(exp_id, user_email, application_url) ).apply_async( link_error=process_error_state.si(exp_id) )

@celery.task
def start_import(exp_id, session_id, user_email, application_url, MAX_BATCH_SIZE = 1000):
    log.debug("starting import")
    
    try:
        upload_helpers.mark_experiment(exp_id, 'loading')
        session = upload.getSessionById(session_id, secure=False)

        accessions, peptides, mod_map, data_runs, errs1 = upload_helpers.parse_datafile(session)
        series_headers = upload_helpers.get_series_headers(session)
        
        log.debug("Getting proteins from NCBI: ")
        prot_map, errs2 = upload_helpers.get_proteins_from_ncbi(accessions, MAX_BATCH_SIZE)
        
        for error in errs1 + errs2:
            experiment.createExperimentError(exp_id, error.row, error.msg)
        

        
        log.debug("Creating subtask tree")
        import_tasks = create_import_tasks(exp_id, prot_map, accessions, peptides, mod_map, series_headers, session.units, data_runs)
        
        
        invoke(import_tasks, exp_id, user_email, application_url)
    except Exception:
        log.debug(traceback.format_exc())
        log.debug("Error in import...")