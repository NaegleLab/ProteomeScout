from celery.task import task
from Bio import Entrez
from ptmscout.config import settings
from Bio import Medline
from ptmscout.utils import protein_utils, uploadutils
from ptmworker import pfam_tools
from geeneus import Proteome
import logging

log = logging.getLogger('ptmscout')

@task(rate_limit='3/s', max_retries=5, default_retry_delay=1)
def get_pubmed_record_by_id(pmid):
    Entrez.email = settings.adminEmail
    handle = Entrez.efetch(db="pubmed",id=pmid,rettype="medline",retmode="text")
    records = Medline.parse(handle)
    
    rec_arr = []
    for r in records:
        rec_arr.append(r)
    
    return rec_arr[0] if len(rec_arr) == 1 else None


def get_qualifier(name, feature_quals):
    qualifier = [ q for q in feature_quals if q['GBQualifier_name'] == name ]
    if len(qualifier) == 1:
        return qualifier[0]['GBQualifier_value']
    return None

def parse_proteinxml_pfam_domains(feature_table):
    parsed_features = []
    for f in feature_table:
        note_val = get_qualifier('note', f['GBFeature_quals'])
        if f['GBFeature_key'] == 'Region' and note_val != None and note_val.find('pfam') > 0:
            domain = pfam_tools.PFamDomain()
            domain.significant = 1
            domain.release = 0
            domain.p_value = -1
            domain.label = get_qualifier('region_name', f['GBFeature_quals'])
            domain.accession = f['GBFeature_intervals'][0]['GBInterval_accession']
            domain.start = int(f['GBFeature_intervals'][0]['GBInterval_from'])
            domain.stop = int(f['GBFeature_intervals'][0]['GBInterval_to'])
            domain.class_ = "Domain"
            
            parsed_features.append(domain)
        
    return parsed_features


def get_gene_name(feature_table):
    for f in feature_table:
        if f['GBFeature_key'] == 'gene':
            gene_val = get_qualifier('gene', f['GBFeature_quals'])
            return gene_val


def parse_sequence_id(seqid):
    result = []
    if seqid.find('sp') == 0:
        _, sp_acc, sp_locus = seqid.split("|")
        i = sp_acc.rfind(".")
        if i > 0:
            sp_acc = sp_acc[0:i]
        
        result.append(('swissprot', sp_acc))
        result.append(('swissprot', sp_locus))
    else:
        seq_type = protein_utils.get_accession_type(seqid)
        result.append((seq_type, seqid))
        
    return result

class EntrezError(Exception):
    def __init__(self, accession):
        self.acc = accession
        
    def __repr__(self):
        return "Unable to fetch protein accession: %s" % (self.acc)

def get_protein_information(pm, acc):
    seq = pm.get_protein_sequence(acc).upper()
    
    if seq == "":
        raise EntrezError(acc)
    
    XML = pm.get_raw_xml(acc)[0]
    
    name = XML['GBSeq_definition']
    gene = get_gene_name(XML['GBSeq_feature-table'])
    taxonomy = XML['GBSeq_taxonomy']
    taxonomy = set([ t.strip().lower() for t in taxonomy.split(';') ])
    
    species = XML['GBSeq_organism']
    
    acc_ids_to_parse = set(XML['GBSeq_other-seqids'])
    acc_ids_to_parse.add(acc)
    acc_ids_to_parse.add( XML['GBSeq_primary-accession'] )
    
    prot_accessions = []
    for seqid in acc_ids_to_parse:
        prot_accessions.extend(parse_sequence_id(seqid))
        
    prot_accessions = set(prot_accessions)        
    prot_domains = parse_proteinxml_pfam_domains(XML['GBSeq_feature-table'])
    
    return name, gene, taxonomy, species, prot_accessions, prot_domains, seq


def load_proteins(accessions, pm):
    log.debug("Querying for proteins: %s", str(accessions))
    pm.batch_get_protein_sequence(accessions)


def get_proteins_from_ncbi(accessions, MAX_BATCH_SIZE):
    pm = Proteome.ProteinManager(settings.adminEmail)
    load_proteins(accessions, pm)

    prot_map = {}
    errors = []
    for acc in accessions:
        try:
            prot_map[acc] = get_protein_information(pm, acc)
        except EntrezError, e:
            errors.append(e)
    
    return prot_map, errors