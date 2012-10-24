import re


def get_accession_type(acc):
    acc_type = None
    
    if(re.search('^gi', acc) != None):
        acc_type = 'gi'
    elif(re.search('^[NXZ]P_\d+', acc) != None):
        acc_type = 'refseq'
    elif(re.search('^[O|P|Q]\d...\d([\.\-]\d+)?$', acc) != None):
        acc_type = 'swissprot'
    elif(re.search('^[A-N|R-Z]\d[A-Z]..\d(\.\d+)?$', acc) != None):
        acc_type = 'swissprot'
    elif(re.search('^[A-Z]{4}_[A-Z]+$', acc) != None):
        acc_type = 'swissprot'
    elif(re.search('^[A-Z]{3}\d{5}$', acc) != None):
        acc_type = 'genbank'
    elif(re.search('^IPI\d+(\.\d+)?$', acc) != None):
        acc_type = 'ipi'
    elif re.search('^ENS', acc) != None:
        acc_type = 'ensembl'
        
    return acc_type
