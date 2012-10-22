import os
import math
from celery.task import task
from ptmscout.config import settings
from celery.canvas import group, chain


@task
def finalize_import():
    pass

@task(rate_limit='3/s') 
def load_proteins(accessions):
    pass

@task
def load_peptide(protein_map, accession, pep_seq):
    pass

@task
def insert_run_data(MSpeptide, series_header, run_name, series):
    pass

@task
def start_import(exp, column_map={}):
    MAX_BATCH_SIZE = 1000
    
    f = open(os.path.join(settings.experiment_data_file_path, exp.datafile), 'rb')
    header = f.readline().strip("\t")
    
    accessions = set()
    peptides = {}
    data_runs = {}
    
    import_tasks = []
    series_headers = [ header[c] for c in column_map['data'] ]
    
    for line in f:
        data = line.split("\t")
        
        acc = data[column_map['accession']]
        accessions.add(acc)
        
        pep = data[column_map['peptide']]
        peps = peptides.get(acc, set())
        peps.add(pep)
        peptides[acc] = peps 
        
        runs = data_runs.get((acc,pep), {})
        run = data[column_map['run']]
        series = []
        
        for c in column_map['data']:
            series.append(float(data[c]))
        
        runs[run] = series
        data_runs[(acc,pep)] = runs
    
    
    num_batches = max([3, int(math.ceil(len(accessions) / float(MAX_BATCH_SIZE)))])
    BATCH_SIZES = [len(accessions) / num_batches] * num_batches
    remainder = len(accessions) - sum(BATCH_SIZES)
    
    for i in xrange(0, remainder):
        BATCH_SIZES[i] += 1
    
    acc_job_args = []
    pep_tasks = []
        
    batch = 0
    for acc in accessions:
        
        acc_job_args.append(acc)
        
        for pep in peptides[acc]:
            run_tasks = []
            
            for name, entry in data_runs[(acc, pep)].items():
                run_tasks.append(insert_run_data.s(series_headers, name, entry))
            
            pep_tasks.append(load_peptide.s(pep))
            pep_tasks.append( chain( load_peptide.s(acc, pep) | group(run_tasks) ) )
        
        if len(acc_job_args) == BATCH_SIZES[batch]:
            batch+=1
            import_tasks.append( chain( load_proteins.s(acc_job_args) | group(pep_tasks) ) )
            acc_job_args = []
            pep_tasks = []

    res = group(import_tasks).apply_async(link=finalize_import.s())
    
    exp.import_process_id = res.id
    exp.status = 'loading'
    exp.saveExperiment()
    
    return res.id