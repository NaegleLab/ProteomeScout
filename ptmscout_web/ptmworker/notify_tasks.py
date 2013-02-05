import celery
import logging
from ptmworker.helpers import upload_helpers
from ptmscout.database import experiment, modifications
from ptmscout.config import strings, settings
from ptmscout.utils import mail

@celery.task
@upload_helpers.transaction_task
def set_loading_status(exp_id, status):
    exp = experiment.getExperimentById(exp_id, check_ready=False, secure=False)
    exp.status = status
    exp.saveExperiment()

@celery.task
@upload_helpers.transaction_task
def set_loading_stage(exp_id, stage, max_value):
    exp = experiment.getExperimentById(exp_id, check_ready=False, secure=False)
    exp.loading_stage = stage
    exp.saveExperiment()

    experiment.setExperimentProgress(exp_id, 0, max_value)

@celery.task
@upload_helpers.transaction_task
def set_progress(exp_id, value, max_value):
    experiment.setExperimentProgress(exp_id, value, max_value)

@celery.task
@upload_helpers.transaction_task
def finalize_experiment_error_state(exc, stack_trace, exp_id, user_email, application_url):
    exp = experiment.getExperimentById(exp_id, check_ready=False, secure=False)

    exp.failure_reason = stack_trace
    exp.status = 'error'
    exp.saveExperiment()

    message = "Exception: " + str(exc)
    subject = strings.experiment_upload_failed_subject
    message = strings.experiment_upload_failed_message % (exp.name, exp.loading_stage, message, application_url)
    
    mail.celery_send_mail([user_email, settings.adminEmail], subject, message)

@celery.task
def finalize_experiment_error_state_callback(uuid, exp_id, user_email, application_url):
    result = celery.result.AsyncResult(uuid)
    exc = result.get(propagate=False)

    finalize_experiment_error_state(exc, str(result.traceback), exp_id, user_email, application_url)

@celery.task
@upload_helpers.transaction_task
def finalize_import(exp_id, user_email, application_url):
    exp = experiment.getExperimentById(exp_id, check_ready=False, secure=False)
    exp.status = 'loaded'
    exp.saveExperiment()

    peptides = modifications.countMeasuredPeptidesForExperiment(exp_id)
    proteins = modifications.countProteinsForExperiment(exp_id)
    exp_errors = experiment.countErrorsForExperiment(exp_id)

    error_log_url = "%s/experiments/%d/errors" % (application_url, exp_id)
    
    subject = strings.experiment_upload_finished_subject
    message = strings.experiment_upload_finished_message % (exp.name, peptides, proteins, exp_errors, error_log_url)
    
    mail.celery_send_mail([user_email], subject, message)



