import celery
import logging
from ptmworker.helpers import upload_helpers
from ptmscout.database import experiment, upload, jobs, annotations,\
    modifications
from ptmscout.utils import uploadutils
from ptmscout.config import strings
from ptmworker import notify_tasks
import traceback
import sys
log = logging.getLogger('ptmscout')

@celery.task
@upload_helpers.notify_job_failed
@upload_helpers.transaction_task
def annotate_experiment(exp_id, job_id):
    exp = experiment.getExperimentById(exp_id, check_ready=False, secure=False)
    exp.modifications = []

    modified_residues, ptms = upload_helpers.summarize_experiment_load(exp.measurements)

    exp.modified_residues = "".join(modified_residues)

    for ptm in ptms:
        exp.modifications.append(ptm)

    exp.saveExperiment()


def parseAnnotationColumns(session_columns):
    ms_col = None
    annotation_cols = []

    for col in session_columns:
        if col.type == 'MS_id':
            ms_col = col
        elif col.type in set(['numeric','nominative','cluster']):
            annotation_cols.append(col)
    
    return ms_col, annotation_cols


def create_annotation(valid_msIds, ms_col, annotation_col, data_rows):
    ms_annotations = []
    errors = []
    for i, row in enumerate(data_rows):
        try:
            ms_id = int(row[ms_col.column_number].strip())
            value = row[annotation_col.column_number].strip()
           
            if ms_id not in valid_msIds:
                raise uploadutils.ParseError(i+1, None, "Invalid MS_id specified: %d" % (ms_id)) 
            
            if value == '':
                value = None
            elif annotation_col.type == 'numeric':
                try:
                    value = float(value)
                except ValueError:
                    value = None
            
            annotation = annotations.Annotation()
            annotation.MS_id = ms_id
            annotation.value = value
            ms_annotations.append(annotation)
            
        except IndexError:
            errors.append(uploadutils.ParseError(i+1, None, strings.experiment_upload_warning_missing_column))
        except ValueError:
            errors.append(uploadutils.ParseError(i+1, None, strings.experiment_upload_warning_columns_values_should_be % 'integer'))
        except uploadutils.ParseError, pe:
            errors.append(pe)
        except Exception, e:
            errors.append(uploadutils.ParseError(i+1, None, "Unexpected error: " + str(e)))
                
    return ms_annotations, errors


def get_valid_msids(job, session):
    exp = experiment.getExperimentById(session.experiment_id, job.user)
    ms_set = set()
    for ms in exp.measurements:
        ms_set.add(ms.id)
    return ms_set

def name_in_use(experiment_id, user, label):
    user_annotations = annotations.getUserAnnotations(experiment_id, user)
    return reduce(bool.__or__, [ label == annotation.name for annotation in user_annotations ], False) 

@celery.task
@upload_helpers.transaction_task
@upload_helpers.notify_job_failed
def start_annotation_import(session_id, job_id):
    notify_tasks.set_job_status.apply_async((job_id, 'running'))
    
    job = jobs.getJobById(job_id)
    session = upload.getSessionById(session_id, job.user)
    
    _, data_rows = uploadutils.load_header_and_data_rows(session.data_file, sys.maxint)
    ms_col, annotation_cols = parseAnnotationColumns(session.columns)
    valid_msIds = get_valid_msids(job, session)
    
    notify_tasks.set_job_stage.apply_async((job_id, 'annotate', len(annotation_cols)))
    
    errors = []
    i = 0
    total_annotations = 0
    for col in annotation_cols:
        if name_in_use(session.experiment_id, job.user, col.label):
            raise Exception("Annotation column label '%s' is already in use" % (col.label))
        annotation_set = annotations.AnnotationSet()
        annotation_set.name = col.label
        annotation_set.owner_id = job.user.id
        annotation_set.experiment_id = session.experiment_id
        annotation_set.type = col.type
        
        ms_annotations, annot_errors = create_annotation(valid_msIds, ms_col, col, data_rows)
        errors.extend( annot_errors )
        
        annotation_set.annotations = ms_annotations
        annotation_set.save()
        
        total_annotations += len(ms_annotations)
        i+=1
        notify_tasks.set_job_progress.apply_async((job_id, i, len(annotation_cols)))
    
    notify_tasks.finalize_annotation_upload_job.apply_async((job_id, total_annotations, errors))
