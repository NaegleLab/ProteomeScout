from celery.task import task
from Bio import Entrez, pairwise2
from ptmscout.config import settings
from Bio import Medline
from ptmworker.helpers import pfam_tools, upload_helpers
from geeneus import Proteome
import logging
import re
from Bio.SubsMat import MatrixInfo
from ptmscout.utils import protein_utils

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


def map_domain_to_sequence(seq1, domain, seq2):
    domain_seq = seq1[domain.start:domain.stop+1]
    
    new_start = seq2.find(domain_seq)
    
    if new_start == -1:
        return None
    
    new_stop = new_start + len(domain_seq) - 1
    
    new_domain = pfam_tools.PFamDomain()
    new_domain.accession = domain.accession
    new_domain.class_ = domain.class_
    new_domain.label = domain.label
    new_domain.start = new_start
    new_domain.stop = new_stop
    new_domain.p_value = domain.p_value
    new_domain.significant = domain.significant
    new_domain.release = domain.release
    
    return new_domain


def parse_pfam_domains(pfam_domains):
    parsed_features = []
    for pfd in pfam_domains:
        domain = pfam_tools.PFamDomain()
        domain.significant = 1
        domain.release = 0
        domain.p_value = -1
        domain.label = pfd['label']
        domain.accession = pfd['accession']
        domain.start = pfd['start']
        domain.stop = pfd['stop']
        domain.class_ = "Domain"
        
        parsed_features.append(domain)
        
    return parsed_features


class EntrezError(Exception):
    def __init__(self):
        pass
        
    def __repr__(self):
        return "Unable to fetch protein accession: %s" % (self.acc)

def get_feature(name, table):
    for entry in table:
        if 'GBFeature_key' in entry and entry['GBFeature_key'] == name:
            return entry

def get_qualifier(name, entry):
    if 'GBFeature_quals' in entry:
        for qual in entry['GBFeature_quals']:
            if qual['GBQualifier_name'] == name:
                return qual

def parse_species_name(name):
    m = re.match(r'(.*) \(.*\)', name)
    if m:
        return m.group(1)
    return name

def parse_organism_host(xml):
    source = get_feature('source', xml['GBSeq_feature-table'])
    if source:
        qual = get_qualifier('host', source)
        if qual:
            return parse_species_name(qual['GBQualifier_value']).strip()

def filter_update_other_accessions(other_accessions):
    type_map = {'ref':'refseq', 'gb':'genbank', 'emb': 'embl', 'uniprotkb/swiss-prot':'swissprot', 'swissprot-locus':'swissprot', 'international protein index':'ipi'}
    filtered_accessions = set()
    for tp, value in other_accessions:
        m = re.match('^([a-zA-Z]+)\|(.*)\|$', value)
        gi_re = re.match('^[0-9]+$', value)
        if m:
            tp = m.group(1).lower()
            val = m.group(2)
            if tp in type_map:
                filtered_accessions.add(( type_map[tp], val ))

        elif tp=='GI' and gi_re:
            tp = tp.lower()
            val = "gi|%s" % (value)
            filtered_accessions.add(( tp, val ))
        elif tp.lower() in type_map:
            filtered_accessions.add(( type_map[tp.lower()], value ))
        else:
            filtered_accessions.add(( tp.lower(), value ))
    return list(filtered_accessions)




def get_protein_information(pm, acc):
    seq = pm.get_protein_sequence(acc).strip().upper()
    
    if seq == "":
        e = EntrezError()
        e.acc = acc
        raise e
    
    name = pm.get_protein_name(acc)
    gene = pm.get_gene_name(acc)
    taxonomy = pm.get_taxonomy(acc)
    species = pm.get_species(acc).strip()
    prot_accessions = pm.get_other_accessions(acc)

    prot_accessions = filter_update_other_accessions(prot_accessions)

    prot_domains = parse_pfam_domains(pm.get_domains(acc))

    prot_mutations = upload_helpers.parse_variants(acc, seq, pm.get_variants(acc))

    host_organism = None
    try:
        host_organism = parse_organism_host(pm.get_raw_xml(acc))
    except:
        pass

    return name, gene, taxonomy, species, host_organism, prot_accessions, prot_domains, prot_mutations, seq

def get_alignment_scores(seq1, seq2):
    matrix = MatrixInfo.blosum62
    gap_open = -10
    gap_extend = -0.5

    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
    top_aln = alns[0]
    aln_seq1, aln_seq2, _score, _begin, _end = top_aln
    
    return aln_seq1.count("-"), aln_seq2.count("-") 

def get_proteins_from_ncbi(accessions):
    pm = Proteome.ProteinManager(settings.adminEmail, uniprotShortcut=False)
 
    query_accessions = accessions
    
    log.info("Querying for proteins: %s", str(query_accessions))
    
    pm.batch_get_protein_sequence(query_accessions)

    prot_map = {}
    errors = []
    for acc in query_accessions:
        try:
            prot_map[acc] = get_protein_information(pm, acc)
        except EntrezError, e:
            errors.append(e)
    
    return prot_map, errors
