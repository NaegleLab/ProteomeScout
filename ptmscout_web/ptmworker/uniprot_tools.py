from ptmscout.config import settings
import urllib2
import MultipartPostHandler as post_handler
from ptmscout.utils import crypto
import os
from Bio import SeqIO
import re

uniprot_url = 'http://www.uniprot.org/uniprot/?'
uniprot_batch_url = 'http://www.uniprot.org/batch/'

def get_isoform_map(accs):
    isoform_map = {}
    new_accs = set()
    
    for acc in accs:
        m = re.match(r"^([A-Za-z0-9]+)\-([0-9]+)$", acc.strip())
        if m:
            root = m.group(1)
            isoform = int(m.group(2))
            
            requested_isoforms = isoform_map.get(root, set())
            requested_isoforms.add(isoform)
            isoform_map[root] = requested_isoforms
            
            new_accs.add(root)
        else:
            new_accs.add(acc)
            
    return list(new_accs), isoform_map
    

def save_accessions(accs):
    tmpfile = os.path.join(settings.ptmscout_scratch_space, crypto.randomString(10))
    with open(tmpfile, 'w') as accfile:
        accfile.write('\n'.join(accs))
    return open(tmpfile)

def parse_name(name):
    m = re.match(r'sp\|([A-Z0-9\-]+)\|(.*)', name)
    
    if m:
        return m.group(1)
    
def parse_isoform_number(acc):
    m = re.search(r"^([A-Za-z0-9]+)\-([0-9]+)$", acc.strip())
    if m:
        return m.group(1), int(m.group(2))
    
    return acc, 0


MAX_RECORD_PER_ISOFORM_QUERY = 20

def __query_for_isoforms(isoform_map, result_map):
    acc_q = '+OR+'.join([ 'accession:%s' % (k) for k in isoform_map.keys() ])
    query = 'query=%s&format=fasta&include=yes' % (acc_q)
    result = urllib2.urlopen(uniprot_url + query)
    parsed_result = SeqIO.parse(result, 'fasta')
    
    for record in parsed_result:
        full_acc = parse_name(record.name)
        root_acc, isoform_number = parse_isoform_number(full_acc)
        
        if isoform_number in isoform_map[root_acc]:
            result_map[full_acc] = (root_acc, isoform_number, record.seq)

def get_some(source_map, quantity):
    result = {}
    map_keys = source_map.keys()
    
    while len(map_keys) > 0 and len(result) < quantity:
        k = map_keys.pop(0)
        result[k] = source_map[k]
        
        del source_map[k]
    
    return result
         

def get_protein_isoforms(isoform_map):
    if settings.DISABLE_UNIPROT_QUERY:
        return []

    map_copy = isoform_map.copy()
    result_map = {}

    while len(map_copy) > 0:
        query_items = get_some(map_copy, MAX_RECORD_PER_ISOFORM_QUERY)
        __query_for_isoforms(query_items, result_map)
        
    return result_map


def get_uniprot_records(accs):
    if settings.DISABLE_UNIPROT_QUERY:
        return []

    root_accs, isoforms = get_isoform_map(accs)
    instream = save_accessions(root_accs)
    
    opener = urllib2.build_opener(post_handler.MultipartPostHandler())
    result = opener.open(uniprot_batch_url, {'file':instream, 'format':'xml', 'include':'yes'})
    
    parsed_result = SeqIO.parse(result, 'uniprot-xml')
    
#    for record in parsed_result:
#        print record