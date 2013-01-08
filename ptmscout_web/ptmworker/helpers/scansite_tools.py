import urllib2
import xml.dom.minidom as xml 
from ptmscout.config import settings
from ptmscout.utils.decorators import rate_limit



scansite_url = "http://scansite3.mit.edu/Scansite3Webservice/proteinScan/proteinName=PTMSCOUT_QUERY/sequence=%s/motifClass=%s/motifNicknames=/stringencyValue=LOW"

class MotifNode(object):
    def __init__(self, siteNode):
        self.name       = self.__get_node_text(siteNode.getElementsByTagName('motifName')[0])
        self.nickname   = self.__get_node_text(siteNode.getElementsByTagName('motifNickName')[0])
        self.score      = float(self.__get_node_text(siteNode.getElementsByTagName('score')[0]))
        self.sequence   = self.__get_node_text(siteNode.getElementsByTagName('siteSequence')[0])

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

@rate_limit(rate=3)
def get_scansite_motif(pep_seq, motif_class):
    if settings.DISABLE_SCANSITE:
        return []
    
    query_url = scansite_url  % (pep_seq, motif_class)
    
    result = urllib2.urlopen(query_url)
    
    parsed_data = ScansiteParser(result.read())
    predicted_sites = [ motif for motif in parsed_data.sites if pep_seq == motif.sequence ]
    
    return predicted_sites
