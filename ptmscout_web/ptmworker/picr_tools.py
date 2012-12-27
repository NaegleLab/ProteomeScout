import xml.dom.minidom as xml
from xml.parsers.expat import ExpatError
import urllib2
from ptmscout.config import settings



class PICRParser(object):

    def __init__(self, picrxml):
        self.dom = xml.parseString(picrxml)
        
        cross_refs = self.dom.getElementsByTagName('identicalCrossReferences')
        self.references = []
        
        for cross_ref in cross_refs:
            self.parseCrossRef(cross_ref)
    
    def __get_node_text(self, node):
        txt = ""
        for c in node.childNodes:
            txt += str(c.nodeValue) 
        return txt
    
    def parseCrossRef(self, cross_ref):
        dbName = self.__get_node_text(cross_ref.getElementsByTagName('databaseName')[0])
        accession = self.__get_node_text(cross_ref.getElementsByTagName('accession')[0])
        version = self.__get_node_text(cross_ref.getElementsByTagName('accessionVersion')[0])
        
        self.references.append((dbName.lower(), accession, version))

picr_accession_query_url = "http://www.ebi.ac.uk/Tools/picr/rest/getUPIForAccession?accession=%s&taxid=%d&database=IPI&database=REFSEQ&database=SWISSPROT"
picr_sequence_query_url  = "http://www.ebi.ac.uk/Tools/picr/rest/getUPIForSequence?sequence=%s&taxid=%d&database=IPI&database=REFSEQ&database=SWISSPROT"

def get_picr(accession, taxon_id):
    if settings.DISABLE_PICR:
        return []
    
    result = urllib2.urlopen(picr_accession_query_url % (accession, taxon_id))
    
    try:
        return PICRParser(result.read()).references
    except ExpatError:
        return []
    
    
