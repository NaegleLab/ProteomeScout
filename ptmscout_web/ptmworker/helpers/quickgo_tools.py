import requests, sys
from requests.exceptions import HTTPError
import logging
import json
from ptmscout.config import settings
from ptmscout.utils.decorators import rate_limit

log = logging.getLogger('ptmscout')



quick_go_term_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%s/ancestors?relations=is_a"
quick_go_annot_url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId="
RETRY_COUNT = 3


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
        self.children = []

def get_GO_term(goId):
    
    requestURL = quick_go_term_url % ( goId )
    tn = TermNode()
    i = 0
    while i < RETRY_COUNT:
        i+=1
        try:
            response = requests.get(requestURL)
            response.raise_for_status()
            result = response.json()['results'][0]

            tn.goId = result['id']
            tn.goName = result['name']
            tn.goFunction = GO_function_map[result['aspect']]
            for child in result['children']:
                if child['relation'] == 'is_a':
                    tn.children.append(child['id'])
            for ancestor in result['ancestors']:
                tn.is_a.append(ancestor)
            break
        except HTTPError as http_err:
            log.error(http_err)
        except Exception as err:
            log.exception(err)
        else:
            log.info("Success! Get Go Term " + goId)
        return -1, tn

def parse_result(response, annotations, gene_symbols):
    
    for result in response.json()['results']:
        proteinId = result['geneProductId'].split(':')[1]
        goId = result['goId']
        dateAdded = int(result['date'])
        geneSymbol = result['symbol']
        evidence = result['goEvidence']
        
        gene_symbols[proteinId] = geneSymbol
        go_terms = annotations.get(proteinId, {})
        
        if evidence != 'IEA' and (goId not in go_terms or go_terms[goId] < dateAdded):
            go_terms[goId] = dateAdded

        annotations[proteinId] = go_terms

def batch_get_GO_annotations(protein_accessions):
    log.info("Getting go annotations for %d accessions", len(protein_accessions))

    annotations = {}
    gene_symbols = {}
    
    for acc in protein_accessions:
        annotations[acc] = {}
        gene_symbols[acc] = ""
    
    if settings.DISABLE_QUICKGO:
        return annotations, gene_symbols

    i = 0
    while i < RETRY_COUNT:
        i+=1
        try:
            requestURL = quick_go_annot_url % (",".join(protein_accessions))
            response = requests.get(requestURL, headers={ "Accept" : "application/json"})
            response.raise_for_status()
            parse_result(response, annotations, gene_symbols)
            break
        except HTTPError as http_err:
            log.error(http_err)
        except Exception as err:
            log.exception(err)
        else:
            log.info("Success! Protein Accessions " + ",".join(protein_accessions))
        

    return annotations, gene_symbols

# import urllib2
# import csv
# import xml.dom.minidom as xml 
# from ptmscout.config import settings
# import logging
# import httplib
# from ptmscout.utils.decorators import rate_limit

# log = logging.getLogger('ptmscout')

# quick_go_term_url  = "http://www.ebi.ac.uk/QuickGO/GTerm?id=%s&format=oboxml"
# quick_go_annot_url = "http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&protein=%s&col=proteinID,proteinSymbol,goID,evidence,date"

# GO_function_map = {'cellular_component':'C',
#                 'molecular_function':'F',
#                 'biological_process':'P'
#                 }



# class TermNode():
#     def __init__(self):
#         self.goId = None
#         self.goName = None
#         self.goFunction = None
#         self.is_a = []

# class OBOXMLParser(object):
#     def __init__(self, annot_stream):
#         self.dom = xml.parseString(annot_stream)
        
#         domnode = self.dom.getElementsByTagName("obo")[0]
#         oboheader = domnode.getElementsByTagName("header")[0]
#         oboterms = domnode.getElementsByTagName("term")
        
#         self.parseHeader(oboheader)
#         self.entries = []
#         for termnode in oboterms:
#             self.entries.append(self.parseNode(termnode))
        
#     def __get_node_text(self, node):
#         txt = ""
#         for c in node.childNodes:
#             txt += str(c.nodeValue) 
#         return txt
    
#     def parseHeader(self, headernode):
#         version_node = headernode.getElementsByTagName('format-version')[0]
#         self.version = self.__get_node_text(version_node).strip()
             
#     def parseNode(self, termnode):
#         tn = TermNode()
        
#         tn.goId = self.__get_node_text(termnode.getElementsByTagName('id')[0]).strip()
#         tn.goName = self.__get_node_text(termnode.getElementsByTagName('name')[0]).strip()
#         tn.goFunction = self.__get_node_text(termnode.getElementsByTagName('namespace')[0]).strip()
        
#         tn.goFunction = GO_function_map[tn.goFunction]
        
#         for is_a_node in termnode.getElementsByTagName('is_a'):
#             tn.is_a.append(self.__get_node_text(is_a_node).strip())

#         return tn
    
# @rate_limit(rate=10)
# def get_GO_term(goId):
#     log.debug("Getting go term: %s", goId)
#     query_url = quick_go_term_url % (goId)
#     annot_stream = urllib2.urlopen(query_url)
    
#     goxml = OBOXMLParser(annot_stream.read())
    
#     return goxml.version, goxml.entries[0]

# RETRY_COUNT = 3

# def parse_result(annot_stream, annotations, gene_symbols):
#     reader = csv.DictReader(annot_stream, delimiter='\t')
    
#     for row in reader:
#         proteinId = row['ID'].strip()
#         goId = row['GO ID'].strip()
#         dateAdded = row['Date'].strip()
#         geneSymbol = row['Symbol'].strip()
#         evidence = row['Evidence'].strip()
        
#         gene_symbols[proteinId] = geneSymbol
#         go_terms = annotations.get(proteinId, {})
        
#         if evidence != 'IEA' and (goId not in go_terms or go_terms[goId] < dateAdded):
#             go_terms[goId] = dateAdded

#         annotations[proteinId] = go_terms

# @rate_limit(rate=3)
# def batch_get_GO_annotations(protein_accessions):
#     log.info("Getting go annotations for %d accessions", len(protein_accessions))

#     annotations = {}
#     gene_symbols = {}
    
#     for acc in protein_accessions:
#         annotations[acc] = {}
#         gene_symbols[acc] = ""
    
#     if settings.DISABLE_QUICKGO:
#         return annotations, gene_symbols

#     i = 0
#     while i < RETRY_COUNT:
#         i+=1
#         try:
#             query_url = quick_go_annot_url % (",".join(protein_accessions))
#             annot_stream = urllib2.urlopen(query_url)
            
#             parse_result(annot_stream, annotations, gene_symbols)
#             break
#         except urllib2.HTTPError, e:
#             log.warning("HTTPError %s when querying QuickGO, retrying (%d / %d)", str(e), i, RETRY_COUNT)
#         except httplib.BadStatusLine:
#             log.warning("Got bad status line from QuickGO, retrying (%d / %d)", i, RETRY_COUNT)
    
#     return annotations, gene_symbols
