import urllib, urllib2
import xml.dom.minidom as xml 
import time
from ptmscout.config import settings
import logging
from ptmscout.database import protein
from xml.parsers.expat import ExpatError
from ptmscout.utils.decorators import rate_limit
import traceback

log = logging.getLogger('ptmscout')
PFAM_DEFAULT_CUTOFF = 0.00001



class PFamDomain(object):
    def __init__(self):
        self.accession = None
        self.label = None
        self.start = None
        self.stop = None
        self.p_value = None
        self.release = None
        self.class_ = None
        self.significant = None
        

class PFamParser(object):
    
    def __init__(self, pfam_xml):
        dom = xml.parseString(pfam_xml).getElementsByTagName('results')[0]
        self.domains = []
        self.parseDom(dom)
        
    def __get_node_text(self, node):
        txt = ""
        for c in node.childNodes:
            txt += str(c.nodeValue) 
        return txt

    def parseDom(self, dom):
        for node in dom.getElementsByTagName('matches'):
            for pnode in node.getElementsByTagName('protein'):
                self.parseProtein(pnode)
    
    def parseProtein(self, dom):
        for dnode in dom.getElementsByTagName('database'):
            self.parseDB(dnode)
    
    def parseDB(self, dom):
        db_release = dom.getAttribute('release')
        
        for node in dom.getElementsByTagName('match'):
            self.parseMatch(node, db_release)
    
    def parseMatch(self, dom, db_release):
        match_id = dom.getAttribute('id')
        match_class = dom.getAttribute('class')
        match_acc =  dom.getAttribute('accession')
        
        for node in dom.getElementsByTagName('location'):
            self.parseLocation(node, match_acc, match_class, match_id, db_release)
        
    def parseLocation(self, dom, accession, class_, label, db_release):
        pep_start = int(dom.getAttribute('start'))
        pep_end = int(dom.getAttribute('end'))
        pvalue = float(dom.getAttribute('evalue'))
        issignificant = int(dom.getAttribute('significant'))
        
        domain = PFamDomain()
        domain.accession = accession
        domain.label = label
        domain.start = pep_start
        domain.stop = pep_end
        domain.p_value = pvalue
        domain.significant = issignificant
        domain.release = db_release
        domain.class_ = class_
        
        self.domains.append(domain)

INTER_QUERY_INTERVAL=1
TIMEOUT=10

class PFamError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __repr__(self):
        return self.msg


def filter_domains(domains):
    domains = [domain for domain in domains if domain.class_=="Domain" and domain.significant==1]
    domains = sorted(domains, key=lambda domain: domain.p_value)
    
    used_sites = set()
    chosen_domains = []
    
    while( len(domains) > 0):
        domain = domains.pop(0)
        range_set = set(range(domain.start, domain.stop+1))
        
        if len(range_set & used_sites) == 0:
            used_sites |= range_set
            chosen_domains.append(domain)
    
    return chosen_domains


PFAM_MIRRORS = ["http://pfam.janelia.org/search/sequence",
                "http://pfam.sanger.ac.uk/search/sequence",
                "http://pfam.sbc.su.se/search/sequence"]

def wait_for_result(jobrequest):
    job_xml = jobrequest.read()
    dom = xml.parseString(job_xml)
    result_url = dom.getElementsByTagName('result_url')[0].childNodes[0].nodeValue
    
    code = 202
    started = time.clock()
    while(code == 202):
        if time.clock() - started > TIMEOUT:
            raise PFamError("Request Timed Out")
        
        resultquery = urllib2.urlopen(result_url)
        code = resultquery.getcode()
        time.sleep(INTER_QUERY_INTERVAL)
 
    return code, resultquery

@rate_limit(rate=3)
def get_computed_pfam_domains(prot_seq, cutoff):
    if settings.DISABLE_PFAM:
        return []
   
    if len(prot_seq) > 10000:
        raise PFamError("Protein Sequence is too long to query")

    args = {'seq':prot_seq, 'output':'xml', 'evalue':cutoff}

    i = 0
    while i < len(PFAM_MIRRORS):
        try:
            jobrequest = urllib2.urlopen(PFAM_MIRRORS[i], urllib.urlencode(args))
            code, resultquery = wait_for_result(jobrequest)

            if code != 200:
                raise PFamError("Got unexpected response code: %d" % (code))

            xmlresult = resultquery.read()
            parsed_pfam = PFamParser(xmlresult)

            return filter_domains(parsed_pfam.domains)
        except ExpatError, e:
            i += 1
            log.warning("PFAM result parsing failed: %s... (%d / %d)", str(e), i, len(PFAM_MIRRORS))
        except urllib2.HTTPError, e:
            i += 1
            log.warning("PFAM query failed, %s... (%d / %d)", str(e), i, len(PFAM_MIRRORS))

    raise PFamError("Unable to query PFam")

def parse_or_query_domains(prot, domains):
    if len(domains) == 0:
        domains = get_computed_pfam_domains(prot.sequence, PFAM_DEFAULT_CUTOFF)
        source = "COMPUTED PFAM"
        params = "pval=%f" % (PFAM_DEFAULT_CUTOFF)
    else:
        source = "PARSED PFAM"
        params = "ALL"

    for domain in domains:
        dbdomain = protein.ProteinDomain()
        dbdomain.p_value = domain.p_value
        dbdomain.start = domain.start
        dbdomain.stop = domain.stop
        dbdomain.source = source
        dbdomain.version = domain.release
        dbdomain.label = domain.label
        dbdomain.params = params
        prot.domains.append(dbdomain)
