import urllib2
import csv
import xml.dom.minidom as xml 
from ptmscout.config import settings

quick_go_term_url  = "http://www.ebi.ac.uk/QuickGO/GTerm?id=%s&format=oboxml"
quick_go_annot_url = "http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&protein=%s&col=proteinID,proteinSymbol,goID,evidence,date"

GO_function_map = {'cellular_component':'C',
                'molecular_function':'F',
                'biological_process':'P'
                }



class TermNode():
    def __init__(self):
        self.goId = None
        self.goName = None
        self.goFunction = None
        self.is_a = []

class OBOXMLParser(object):
    def __init__(self, annot_stream):
        self.dom = xml.parseString(annot_stream)
        
        domnode = self.dom.getElementsByTagName("obo")[0]
        oboheader = domnode.getElementsByTagName("header")[0]
        oboterms = domnode.getElementsByTagName("term")
        
        self.parseHeader(oboheader)
        self.entries = []
        for termnode in oboterms:
            self.entries.append(self.parseNode(termnode))
        
    def __get_node_text(self, node):
        txt = ""
        for c in node.childNodes:
            txt += str(c.nodeValue) 
        return txt
    
    def parseHeader(self, headernode):
        version_node = headernode.getElementsByTagName('format-version')[0]
        self.version = self.__get_node_text(version_node).strip()
             
    def parseNode(self, termnode):
        tn = TermNode()
        
        tn.goId = self.__get_node_text(termnode.getElementsByTagName('id')[0]).strip()
        tn.goName = self.__get_node_text(termnode.getElementsByTagName('name')[0]).strip()
        tn.goFunction = self.__get_node_text(termnode.getElementsByTagName('namespace')[0]).strip()
        
        tn.goFunction = GO_function_map[tn.goFunction]
        
        for is_a_node in termnode.getElementsByTagName('is_a'):
            tn.is_a.append(self.__get_node_text(is_a_node).strip())

        return tn
    
def get_GO_term(goId):
    query_url = quick_go_term_url % (goId)
    annot_stream = urllib2.urlopen(query_url)
    
    goxml = OBOXMLParser(annot_stream.read())
    
    return goxml.version, goxml.entries[0]

def batch_get_GO_annotations(protein_accessions):
    annotations = {}
    gene_symbols = {}
    
    if settings.DISABLE_QUICKGO:
        for acc in protein_accessions:
            annotations[acc] = {}
            gene_symbols[acc] = ""
        return annotations, gene_symbols
            
            
    
    query_url = quick_go_annot_url % (",".join(protein_accessions))
    annot_stream = urllib2.urlopen(query_url)
    
    reader = csv.DictReader(annot_stream, delimiter='\t')
    
    for row in reader:
        proteinId = row['ID'].strip()
        goId = row['GO ID'].strip()
        dateAdded = row['Date'].strip()
        geneSymbol = row['Symbol'].strip()
        evidence = row['Evidence'].strip()
        
        gene_symbols[proteinId] = geneSymbol
        go_terms = annotations.get(proteinId, {})
        
        if evidence != 'IEA' and (goId not in go_terms or go_terms[goId] < dateAdded):
            go_terms[goId] = dateAdded

        annotations[proteinId] = go_terms
    
    return annotations, gene_symbols