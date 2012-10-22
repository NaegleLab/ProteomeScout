from celery.task import task
from ptmscout.config import settings
import os
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
    
    acc_job_args = []
    pep_tasks = []
        
    for acc in accessions:
        
        acc_job_args.append(acc)
        
        for pep in peptides[acc]:
            run_tasks = []
            
            for name, entry in data_runs[(acc, pep)].items():
                run_tasks.append(insert_run_data.s(series_headers, name, entry))
            
            pep_tasks.append(load_peptide.s(pep))
            pep_tasks.append( chain( load_peptide.s(acc, pep) | group(run_tasks) ) )
        
        if len(acc_job_args) == MAX_BATCH_SIZE:
            import_tasks.append( chain( load_proteins.s(acc_job_args) | group(pep_tasks) ) )
            acc_job_args = []
            pep_tasks = []

    if len(acc_job_args) > 0:
        import_tasks.append( chain( load_proteins.s(acc_job_args) | group(pep_tasks) ) )
            
    res = group(import_tasks).apply_async(link=finalize_import.s())
    
    exp.import_process_id = res.id
    exp.status = 'loading'
    exp.saveExperiment()
    
    return res.id