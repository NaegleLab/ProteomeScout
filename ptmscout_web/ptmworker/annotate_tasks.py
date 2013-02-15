import celery
import logging
from ptmworker.helpers import upload_helpers
from ptmscout.database import experiment
import transaction
import traceback
log = logging.getLogger('ptmscout')

@celery.task
@upload_helpers.transaction_task
def annotate_experiment(exp_id):
    exp = experiment.getExperimentById(exp_id, check_ready=False, secure=False)
    exp.modifications = []

    modified_residues, ptms = upload_helpers.summarize_experiment_load(exp.measurements)

    exp.modified_residues = "".join(modified_residues)

    for ptm in ptms:
        exp.modifications.append(ptm)

    exp.saveExperiment()
