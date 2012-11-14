import urllib2
import xml.dom.minidom as xml 

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

class ScansiteParser():
    def __init__(self, scansite_stream):
        self.dom = xml.parseString(scansite_stream)
        self.sites = []
        
        for siteNode in self.dom.getElementsByTagName('predictedSite'):
            mn = MotifNode(siteNode)
            self.sites.append(mn)

def get_scansite_motif(pep_seq, motif_class):
    query_url = scansite_url  % (pep_seq, motif_class)
    
    result = urllib2.urlopen(query_url)
    
    parsed_data = ScansiteParser(result.read())
    predicted_sites = [ motif for motif in parsed_data.sites if pep_seq == motif.sequence ]
    
    return predicted_sites