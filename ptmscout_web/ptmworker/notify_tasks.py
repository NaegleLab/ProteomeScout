import celery
from ptmworker.helpers import upload_helpers
from ptmscout.database import experiment, modifications, jobs
from ptmscout.config import strings, settings
from ptmscout.utils import mail
import traceback

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

@celery.task
@upload_helpers.transaction_task
def finalize_annotation_upload_job(job_id, total, errors):
    job = jobs.getJobById(job_id)
    job.finish()
    job.save()
    
    subject = strings.annotation_upload_finished_subject
    message = strings.annotation_upload_finished_message % (job.name, total, len(errors), job.result_url)
    
    for err in errors:
        message += "%s\n" % ( err.message ) 
    
    mail.celery_send_mail([job.user.email], subject, message)
    

@celery.task
@upload_helpers.transaction_task
def notify_job_failed(job_id, exc, stack_trace):
    job = jobs.getJobById(job_id)
    job.fail(stack_trace)
    job.save()
    
    subject = strings.job_failed_subject
    message = strings.job_failed_message % (job.name, job.stage, "Exception: " + str(exc))
    
    mail.celery_send_mail([job.user.email, settings.adminEmail], subject, message)
    
    

@celery.task
@upload_helpers.transaction_task
def set_job_status(job_id, status):
    job = jobs.getJobById(job_id)
    job.status = status
    job.save()

@celery.task
@upload_helpers.transaction_task
def set_job_stage(job_id, stage, max_value):
    job = jobs.getJobById(job_id)
    job.stage = stage
    job.max_progress = max_value
    job.save()

@celery.task
@upload_helpers.transaction_task
def set_job_progress(job_id, value, max_value):
    job = jobs.getJobById(job_id)
    job.progress = value
    job.max_progress = max_value
    job.save()

def handle_errors(exp_id_arg):
    def decorate(fn):
        def ttask(*args):
            try:
                return fn(*args)
            except Exception, exc:
                exp_id = args[exp_id_arg]
                finalize_experiment_error_state.apply_async((exc, traceback.format_exc(), exp_id))
                raise

        ttask.__name__ = fn.__name__
        return ttask
    return decorate


