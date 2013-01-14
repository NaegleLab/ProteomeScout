from ptmscout.config import settings
from ptmscout.utils.decorators import rate_limit
import urllib2
import httplib
import MultipartPostHandler as post_handler
from ptmscout.utils import crypto
import os
from Bio import SeqIO
import re
import traceback
import logging
import ptmscout.database.taxonomies
import xml.dom.minidom as xml
from xml.parsers.expat import ExpatError
import mmap

log = logging.getLogger('ptmscout')

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

def parse_description(desc):
    m = re.search(r"[Ii]soform \S+", desc)
    if m:
        return m.group(0)

MAX_RECORD_PER_ISOFORM_QUERY = 20

def __query_for_isoforms(root_accessions, result_map):
    acc_q = '+OR+'.join([ 'accession:%s' % (k) for k in root_accessions ])
    query = uniprot_url + 'query=%s&format=fasta&include=yes' % (acc_q)
    result = urllib2.urlopen(query)
    parsed_result = SeqIO.parse(result, 'fasta')

    for record in parsed_result:
        name = parse_description(record.description)
        full_acc = parse_name(record.name)
        root_acc, isoform_number = parse_isoform_number(full_acc)

        if isoform_number > 0:
            result_map[full_acc] = (root_acc, name, isoform_number, record.seq)


def get_protein_isoforms(root_accessions):
    if settings.DISABLE_UNIPROT_QUERY:
        return {}

    result_map = {}

    while len(root_accessions) > 0:
        end = min(MAX_RECORD_PER_ISOFORM_QUERY, len(root_accessions))
        query_items = root_accessions[0:end]
        root_accessions = root_accessions[end:]

        __query_for_isoforms(query_items, result_map)

    return result_map

def get_scientific_name(name):
    m = re.match(r"^(.*) \(.*\)$", name)
    if m != None:
        return m.group(1).strip()
    return name

def parse_xml(xml):
    name = xml.description
    gene = None

    if 'gene_name_primary' in xml.annotations:
        gene = xml.annotations['gene_name_primary']

    taxons = xml.annotations['taxonomy']
    species = get_scientific_name(xml.annotations['organism'])

    other_accessions = [('swissprot', xml.id), ('swissprot', xml.name)]
    if 'gene_name_synonym' in xml.annotations:
        for gene_name in xml.annotations['gene_name_synonym']:
            other_accessions.append(('gene_synonym', gene_name))

    domains = []
    seq = xml.seq

    host_organism = None
    if 'organism_host' in xml.annotations:
        host_organism = xml.annotations['organism_host'][0]

    return xml.id, (name, gene, taxons, species, host_organism, other_accessions, domains, seq)


def handle_result(result):
    try:
        str_result = result.read()
    except httplib.IncompleteRead, e:
        str_result = e.partial

    handle = mmap.mmap(-1, len(str_result))
    handle.write(str_result)
    handle.seek(0)

    parsed_result = SeqIO.parse(handle, 'uniprot-xml')

    result_map = {}
    i = 0
    for xml in parsed_result:
        try:
            acc, prot_info = parse_xml(xml)
            result_map[acc] = prot_info
        except:
            pass
        i+=1

    return result_map


def map_isoform_results(result_map, isoform_map):
    isoforms_by_root = {}

    for iso_acc in isoform_map:
        root_acc, isoform_name, iso_number, iso_seq = isoform_map[iso_acc]

        identified_isoforms = isoforms_by_root.get(root_acc, set())
        isoforms_by_root[root_acc] = identified_isoforms | set([iso_number])

        name, gene, taxons, species, host_organism,  _o, _d, _s = result_map[root_acc]
      
        isoform_fullname = "%s (%s)" % (name, isoform_name)
        isoform_accs = [('swissprot', iso_acc)]
        result_map[iso_acc] = (isoform_fullname, gene, taxons, species, host_organism, isoform_accs, [], iso_seq)

    for root in isoforms_by_root:
        isos = isoforms_by_root[root]
        expected_isos = set( range(1, len(isos)+2) )

        missing_isos = expected_isos - isos
        if len(missing_isos) == 1:
            canonical_isoform = missing_isos.pop()

            new_isoform = "%s-%d" % (root, canonical_isoform)

            isoform_fullname, gene, taxons, species, host_organism, isoform_accs, domains, iso_seq = result_map[root]
            isoform_accs.append(('swissprot', new_isoform))
            
            result_map[root] = isoform_fullname, gene, taxons, species, host_organism, isoform_accs, domains, iso_seq
            result_map[new_isoform] = isoform_fullname, gene, taxons, species, host_organism, isoform_accs, domains, iso_seq

def map_combine(r1, r2):
    return dict(r1.items() + r2.items())

MAX_RETRIES = 3

@rate_limit(rate=3)
def get_uniprot_records(accs):
    log.debug("Query: %s", str(accs))
    if settings.DISABLE_UNIPROT_QUERY:
        return []

    root_accs, isoforms = get_isoform_map(accs)

    i = 0
    while i < MAX_RETRIES:
        try:
            try:
                instream = save_accessions(root_accs)
                opener = urllib2.build_opener(post_handler.MultipartPostHandler())
                result = opener.open(uniprot_batch_url, {'file':instream, 'format':'xml', 'include':'yes'})

                result_map = handle_result(result)
                isoform_map = get_protein_isoforms(isoforms.keys())

                map_isoform_results(result_map, isoform_map)
                return result_map
            except urllib2.HTTPError:
                log.debug( "Failed query..." )
                if len(accs) == 1:
                    return {}
                else:
                    bisect = len(accs) / 2
                    r1 = get_uniprot_records(accs[:bisect])
                    r2 = get_uniprot_records(accs[bisect:])

                    return map_combine(r1, r2)
        except:
            i+=1
            log.info("Uniprot query failed (retry %d / %d)", i, MAX_RETRIES)

    return {}
