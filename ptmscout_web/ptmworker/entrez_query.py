from celery.task import task
from Bio import Entrez, SeqIO
from ptmscout.config import settings
from Bio import Medline

@task(rate_limit='3/s', max_retries=5, default_retry_delay=1)
def get_protein_records_by_accession(accessions):
    Entrez.email = settings.adminEmail
    search_handle = Entrez.esearch(db="protein", term=",".join(accessions), usehistory="y")
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    gi_list = search_results["IdList"]
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"] 
    
    fetch_handle = Entrez.efetch(db="protein", rettype='gb', retmode='text', id=gi_list, webenv=webenv, query_key=query_key)
    record_iterator = SeqIO.parse(fetch_handle, "genbank")
    
    record_map = {}
    try:
        for r in record_iterator:
            for acc in r.annotations['accessions']:
                records = record_map.get(acc, set())
                records.add(r)
                record_map[acc] = records
    except IndexError:
        pass
    
    return record_map
    
@task(rate_limit='3/s', max_retries=5, default_retry_delay=1)
def get_pubmed_record_by_id(pmid):
    Entrez.email = settings.adminEmail
    handle = Entrez.efetch(db="pubmed",id=pmid,rettype="medline",retmode="text")
    records = Medline.parse(handle)
    
    rec_arr = []
    for r in records:
        rec_arr.append(r)
    
    return rec_arr[0] if len(rec_arr) == 1 else None