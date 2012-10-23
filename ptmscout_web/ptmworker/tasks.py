import os
import math
from celery.task import task
from ptmscout.config import settings
from celery.canvas import group, chain
from ptmscout.database import experiment, modifications


@task
def finalize_import(exp):
    from ptmscout.database import DBSession
    exp.status = 'loaded'
    exp.saveExperiment()
    DBSession.flush()
    
    return 0

@task(rate_limit='3/s') 
def load_proteins(accessions):
    pass

@task
def load_peptide(protein_map, exp_id, accession, pep_seq):
    prot = protein_map[accession]
    prot_seq = prot.sequence
    upper_case = pep_seq.upper()
    
    index = prot_seq.find(upper_case)
    
    if index == -1:
        raise Exception("Peptide sequence not found in protein")
    
    mod_sites = []
    for j in xrange(0, len(pep_seq)):
        if pep_seq[j] != upper_case[j]:
            mod_sites.append(j)
    
    MS_pep = modifications.MeasuredPeptide()
    MS_pep.experiment_id = exp_id
    MS_pep.phosphopep = pep_seq
    MS_pep.protein_id = prot.id
    
    
    for i in mod_sites:
        pep_site = i + index
        pep_tryps   = upper_case[:i] + pep_seq[i] + upper_case[i+1:]
        pep_aligned = prot_seq[pep_site-7:pep_site] + pep_seq[i] + prot_seq[pep_site+1:pep_site+8]
        pep_type = upper_case[i]
        
        pep = modifications.Phosphopep()
        
        pep.pep_aligned = pep_aligned
        pep.pep_tryps = pep_tryps
        pep.site_pos = pep_site + 1
        pep.site_type = pep_type
        pep.protein_id = prot.id
        pep.pfam_site = "~~~"
        
        MS_pep.phosphopeps.append(pep)
        
    MS_pep.save()
    
    return MS_pep
    

@task
def insert_run_data(MSpeptide, series_header, run_name, series):
    for i in xrange(0, len(series_header)):
        data = experiment.ExperimentData()
        
        _, t, x = series_header[i].split(":")
        y = float(series[i])
        
        data.run = run_name
        data.priority = i + 1
        data.type = t
        data.label = x
        data.value = y
        data.MS_id = MSpeptide.id
        
        data.save()
        
@task
def start_import(exp, column_map={}, MAX_BATCH_SIZE = 1000):
    f = open(os.path.join(settings.experiment_data_file_path, exp.datafile), 'rb')
    header = f.readline().strip("\t")
    
    accessions = set()
    peptides = {}
    data_runs = {}
    
    import_tasks = []
    series_headers = [ header[c] for c in column_map['data'] ]
    
    for line in f:
        data = line.split("\t")
        
        acc = data[column_map['accession']].strip()
        accessions.add(acc)
        
        pep = data[column_map['peptide']].strip()
        peps = peptides.get(acc, set())
        peps.add(pep)
        peptides[acc] = peps 
        
        runs = data_runs.get((acc,pep), {})
        run = data[column_map['run']].strip()
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
            pep_tasks.append( chain( load_peptide.s(exp.id, acc, pep) | group(run_tasks) ) )
        
        if len(acc_job_args) == BATCH_SIZES[batch]:
            batch+=1
            import_tasks.append( chain( load_proteins.s(acc_job_args) | group(pep_tasks) ) )
            acc_job_args = []
            pep_tasks = []

    res = group(import_tasks).apply_async(link=finalize_import.s(exp))
    
    exp.import_process_id = res.id
    exp.status = 'loading'
    exp.saveExperiment()
    
    return res.id