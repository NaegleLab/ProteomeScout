import urllib2
import httplib
import xml.dom.minidom as xml
from ptmscout.config import settings
from ptmscout.utils.decorators import rate_limit
import logging
from collections import defaultdict
from requests.exceptions import HTTPError
import requests
from xml.etree import ElementTree as ET

log = logging.getLogger('ptmscout')
# scansite_url = "http://scansite3.mit.edu/ws/proteinScan/proteinName=PTMSCOUT_QUERY/sequence=%s/motifClass=%s/motifNicknames=/stringencyValue=LOW"
#scansite_url = "http://scansite3.mit.edu/Scansite3Webservice/proteinScan/proteinName=PTMSCOUT_QUERY/sequence=%s/motifClass=%s/motifNicknames=/stringencyValue=LOW"
scansite_url = "https://scansite4.mit.edu/webservice/proteinscan/identifier=PTMSCOUT_QUERY/sequence=%s/motifclass=%s/stringency=Low"

class MotifNode():
    def __init__(self, siteNode):
        self.name       = siteNode.find('motifName').text
        self.nickname   = siteNode.find('motifShortName').text
        self.score      = siteNode.find('score').text
        self.percentile = siteNode.find('percentile').text
        self.sequence   = siteNode.find('siteSequence').text
        self.site       = siteNode.find('site').text


    def parse_source(self):
        lower_name = self.name.lower()
        lower_nick = self.nickname.lower()
        is_kinase = lower_name.find("kinase") > -1 or lower_nick.find("_kin") > -1
        if is_kinase:
            return 'scansite_kinase'
        return 'scansite_bind'

class ScansiteParser():
    def __init__(self, scansite_stream):
        self.root = ET.fromstring(scansite_stream)
        self.sites = []

        for siteNode in self.root.iter('predictedSite'):
            mn = MotifNode(siteNode)
            self.sites.append(mn) 



class ScansiteError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        return self.msg

RETRY_COUNT = 3

@rate_limit(rate=settings.SCANSITE_RATE_LIMIT)
def get_scansite_motif(pep_seq, motif_class, filter_exact=True):
    if settings.DISABLE_SCANSITE:
        return []
    
    i = 0
    while i < RETRY_COUNT:
        i+=1
        try:
            response = requests.get(scansite_url % (pep_seq, motif_class))
            response.raise_for_status()
            parsed_data = ScansiteParser(response.content)  
            if filter_exact:
                predicted_sites = [ motif for motif in parsed_data.sites if pep_seq == motif.sequence ]
            else:
                predicted_sites = parsed_data.sites

            return predicted_sites
        
        except HTTPError as http_err:
            log.error(http_err)
        except Exception as err:
            log.exception(err)
        
    raise ScansiteError("Failed query: %s" % (query_url))


def chop_overlapping(seq, size, overlap):
    query_sequences = []
    offset = 0
    while( len(seq) > 0 ):
        query_sequences.append( (offset, seq[0:size]) )

        offset += size - overlap
        seq = seq[(size - overlap):]

    return query_sequences

def get_scansite_protein_motifs(prot_seq, motif_class):
    scansite_results = defaultdict(dict)
    for offset, pep_seq in chop_overlapping(prot_seq, 1000, 50):
        predictions = get_scansite_motif(pep_seq, motif_class, filter_exact=False)

        for p in predictions:
            site_type = p.site[0]
            site_pos = int(p.site[1:]) + offset
            p.site = "%s%d" % (site_type, site_pos)
            scansite_results[site_pos][p.nickname] = p

    return [ scansite_results[site][name] for site in sorted(scansite_results.keys()) for name in scansite_results[site] ]
