from celery.task import task
from ptmscout.config import settings
import os

@task 
def load_protein(accession):
    pass

@task
def load_peptide(protein, pep_seq):
    pass

@task
def insert_run_data(protein_id, peptide_id, run_data):
    pass

@task
def start_import(exp, column_map={}):
    open(os.path.join(settings.experiment_data_file_path, exp.dataset), 'rb')