from celery.task import task
from Bio import Entrez, SeqIO
from ptmscout.config import settings
from Bio import Medline
    
@task(rate_limit='3/s', max_retries=5, default_retry_delay=1)
def get_pubmed_record_by_id(pmid):
    Entrez.email = settings.adminEmail
    handle = Entrez.efetch(db="pubmed",id=pmid,rettype="medline",retmode="text")
    records = Medline.parse(handle)
    
    rec_arr = []
    for r in records:
        rec_arr.append(r)
    
    return rec_arr[0] if len(rec_arr) == 1 else None