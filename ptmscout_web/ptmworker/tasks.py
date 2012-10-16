from celery.task import task
from ptmscout.config import settings
import os
from celery.canvas import group, chain, chord

@task
def finalize_import():
    pass

@task 
def load_protein(accession):
    pass

@task
def load_peptide(protein, pep_seq):
    pass

@task
def insert_run_data(MSpeptide, series_header, run_name, series):
    pass

@task
def start_import(exp, column_map={}):
    
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
        
    for acc in accessions:
        pep_tasks = []
        
        for pep in peptides[acc]:
            run_tasks = []
            
            for name, entry in data_runs[(acc, pep)].items():
                
                run_tasks.append(insert_run_data.s(series_headers, name, entry))
            
            pep_tasks.append(load_peptide.s(pep))
            pep_tasks.append( chain( load_peptide.s(pep) | group(run_tasks) ) )
        
        import_tasks.append( chain( load_protein.s(acc) | group(pep_tasks) ) )
     
    res = group(import_tasks).apply_async(link=finalize_import.s())
    
    exp.import_process_id = res.id
    exp.saveExperiment()
    
    return res.id