from celery.task import task
from Bio import Entrez, pairwise2
from ptmscout.config import settings
from Bio import Medline
from ptmworker import pfam_tools
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
            return parse_species_name(qual['GBQualifier_value'])


def get_protein_information(pm, acc):
    seq = pm.get_protein_sequence(acc).upper()
    
    if seq == "":
        e = EntrezError()
        e.acc = acc
        raise e
    
    name = pm.get_protein_name(acc)
    gene = pm.get_gene_name(acc)
    taxonomy = pm.get_taxonomy(acc)
    species = pm.get_species(acc)
    prot_accessions = pm.get_other_accessions(acc)
    prot_domains = parse_pfam_domains(pm.get_domains(acc))

    host_organism = parse_organism_host(pm.get_raw_xml(acc))
    prot_isoforms = []#pm.get_isoforms(acc)
    
    for i in prot_isoforms:
        prot_isoforms[i] = prot_isoforms[i].upper()
    
    return name, gene, taxonomy, species, host_organism, prot_accessions, prot_domains, seq, prot_isoforms


def parse_isoform_number(acc):
    m = re.search(r"^(.*)\-([0-9]+)$", acc.strip())
    if m:
        return m.group(1), m.group(2)
    
    return acc, None

def get_isoform_map(accs):
    isoform_map = {}
    new_accs = set()
    
    for acc in accs:
        root, isoform = parse_isoform_number(acc)
        
        if isoform:
            isoform_map[acc] = root
            new_accs.add(root)
        else:
            new_accs.add(acc)
        
    return list(new_accs), isoform_map

def get_alignment_scores(seq1, seq2):
    matrix = MatrixInfo.blosum62
    gap_open = -10
    gap_extend = -0.5

    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
    top_aln = alns[0]
    aln_seq1, aln_seq2, _score, _begin, _end = top_aln
    
    return aln_seq1.count("-"), aln_seq2.count("-") 


def create_isoform(root_acc, isoform_accession, isoform_id, isoform_seq, name, gene, taxonomy, species, host_organism, prot_domains, seq):
    isoform_domains = []
    isoform_name = "%s isoform %s" % (name, isoform_id)
    inserted, deleted = get_alignment_scores(seq, isoform_seq)
    
    if inserted < settings.isoform_sequence_diff_pfam_threshold and \
        deleted < settings.isoform_sequence_diff_pfam_threshold:
        
        for domain in prot_domains:
            ndomain = map_domain_to_sequence(seq, domain, isoform_seq)
            
            if ndomain == None:
                isoform_domains = []
                break
            
            isoform_domains.append(ndomain)
        
    acc_type = protein_utils.get_accession_type(root_acc)
    return (isoform_name, gene, taxonomy, species, host_organism, [(acc_type, isoform_accession)], isoform_domains, isoform_seq)


def get_proteins_from_ncbi(accessions):
    pm = Proteome.ProteinManager(settings.adminEmail)
 
    query_accessions = accessions
#    query_accessions, _ = get_isoform_map(accessions)
    
    log.info("Querying for proteins: %s", str(query_accessions))
    
    pm.batch_get_protein_sequence(query_accessions)

    prot_map = {}
    errors = []
    for acc in query_accessions:
        try:
            name, gene, taxonomy, species, host_organism, prot_accessions, prot_domains, seq, prot_isoforms = get_protein_information(pm, acc)
            prot_map[acc] = (name, gene, taxonomy, species, host_organism, prot_accessions, prot_domains, seq)
            
            for isoform_id in prot_isoforms:
                isoform_seq = prot_isoforms[isoform_id]
                isoform_accession = "%s-%s" % (acc, isoform_id)
                prot_map[isoform_accession] = create_isoform(acc, isoform_accession, isoform_id, isoform_seq, name, gene, taxonomy, species, host_organism, prot_domains, seq)

        except EntrezError, e:
            errors.append(e)
    
    return prot_map, errors
