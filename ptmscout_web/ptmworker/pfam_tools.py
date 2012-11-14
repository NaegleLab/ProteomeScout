import urllib, urllib2
import xml.dom.minidom as xml 
import time

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
        domain.end = pep_end
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
        range_set = set(range(domain.start, domain.end+1))
        
        if len(range_set & used_sites) == 0:
            used_sites |= range_set
            chosen_domains.append(domain)
    
    return chosen_domains


def get_computed_pfam_domains(prot_seq, cutoff):
    args = {'seq':prot_seq, 'output':'xml', 'evalue':cutoff}
    
    jobrequest = urllib2.urlopen("http://pfam.sanger.ac.uk/search/sequence", urllib.urlencode(args))
    
    dom = xml.parseString(jobrequest.read())
    result_url = dom.getElementsByTagName('result_url')[0].childNodes[0].nodeValue
    
    code = 202
    started = time.clock()
    while(code == 202):
        if time.clock() - started > TIMEOUT:
            raise PFamError("Request Timed Out")
        
        resultquery = urllib2.urlopen(result_url)
        code = resultquery.getcode()
        time.sleep(INTER_QUERY_INTERVAL)
    
    if code != 200:
        raise PFamError("Got unexpected response code: %d" % (code))
    
    parsed_pfam = PFamParser(resultquery.read())
    
    return filter_domains(parsed_pfam.domains)

