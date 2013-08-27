import urllib2
import httplib
import xml.dom.minidom as xml
from ptmscout.config import settings
from ptmscout.utils.decorators import rate_limit
import logging
from collections import defaultdict

log = logging.getLogger('ptmscout')
scansite_url = "http://scansite3.mit.edu/ws/proteinScan/proteinName=PTMSCOUT_QUERY/sequence=%s/motifClass=%s/motifNicknames=/stringencyValue=LOW"
#scansite_url = "http://scansite3.mit.edu/Scansite3Webservice/proteinScan/proteinName=PTMSCOUT_QUERY/sequence=%s/motifClass=%s/motifNicknames=/stringencyValue=LOW"

class MotifNode(object):
    def __init__(self, siteNode):
        self.name       = self.__get_node_text(siteNode.getElementsByTagName('motifName')[0])
        self.nickname   = self.__get_node_text(siteNode.getElementsByTagName('motifNickName')[0])
        self.score      = float(self.__get_node_text(siteNode.getElementsByTagName('score')[0]))
        self.percentile = 100.0 * float(self.__get_node_text(siteNode.getElementsByTagName('percentile')[0]))
        self.sequence   = self.__get_node_text(siteNode.getElementsByTagName('siteSequence')[0])
        self.site       = self.__get_node_text(siteNode.getElementsByTagName('site')[0])

    def __get_node_text(self, node):
        txt = ""
        for c in node.childNodes:
            txt += str(c.nodeValue) 
        return txt

    def parse_source(self):
        lower_name = self.name.lower()
        lower_nick = self.nickname.lower()
        is_kinase = lower_name.find("kinase") > -1 or lower_nick.find("_kin") > -1
        if is_kinase:
            return 'scansite_kinase'
        return 'scansite_bind'

class ScansiteParser():
    def __init__(self, scansite_stream):
        self.dom = xml.parseString(scansite_stream)
        self.sites = []
        
        for siteNode in self.dom.getElementsByTagName('predictedSite'):
            mn = MotifNode(siteNode)
            self.sites.append(mn)

class ScansiteError(Exception):
    pass

RETRY_COUNT = 3

@rate_limit(rate=settings.SCANSITE_RATE_LIMIT)
def get_scansite_motif(pep_seq, motif_class, filter_exact=True):
    if settings.DISABLE_SCANSITE:
        return []

    i = 0
    while i < RETRY_COUNT:
        i+=1
        try:
            query_url = scansite_url  % (pep_seq.strip(), motif_class)
            result = urllib2.urlopen(query_url)

            parsed_data = ScansiteParser(result.read())
            if filter_exact:
                predicted_sites = [ motif for motif in parsed_data.sites if pep_seq == motif.sequence ]
            else:
                predicted_sites = parsed_data.sites

            return predicted_sites

        except urllib2.HTTPError, e:
            log.warning("HTTPError %s when querying Scansite3, retrying (%d / %d)", str(e), i, RETRY_COUNT)
        except httplib.BadStatusLine:
            log.warning("Got bad status line from Scansite3, retrying (%d / %d)", i, RETRY_COUNT)

    raise ScansiteError()

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
