

import pandas as pd
import requests
import re
from io import StringIO
from ptmscout.config import settings

def pull_xref(accession, from_database = 'P_REFSEQ_AC', to_database = 'ACC', xref_db_name = 'swissprot'):
    """
    Pulls Uniprot IDs given the provided RefSeq protein id(s) using Uniprot API
    
    Returns
    --------
    xrefs : list of lists
        database
        xref accession
    """
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
        'from' : from_database,
        'to' : to_database,
        'format' : 'tab',
        'query' : accession
    }
    response = requests.get(url, params)
    if response.ok and response.text:
        xrefs = pd.read_csv(StringIO(response.text), sep = '\t')
        xrefs['From'] = xref_db_name
        return xrefs.values.tolist()
    return []


def get_picr(accession, taxon_id = None, database = None):
    """
    Cross reference accession against Uniprot and RefSeq protein ids
    If database is provided pull using said database
    If no database is provided check against Uniprot and RefSeq databases
    
    Parameters
    -----------
    accession : str
        protein id
    taxon_id : str
        taxon (not currently used)
    database : str
        what database accession maps to
        see https://www.uniprot.org/help/api_idmapping to find database mapping
    
    Returns
    -------
    xrefs : list of lists
        database
        xref accession
    """
    if settings.DISABLE_PICR:
        return []
    
    if database:
        full_xrefs = []
        full_xrefs.extend(pull_xref(accession, from_database = database, to_database = 'ACC', xref_db_name = 'swissprot'))
        full_xrefs.extend(pull_xref(accession, from_database = database, to_database = 'P_REFSEQ_AC', xref_db_name = 'refseq'))
        return full_xrefs
    else:
        xrefs = pull_xref(accession, from_database = 'ACC', to_database = 'P_REFSEQ_AC', xref_db_name = 'refseq')
        if xrefs:
            return xrefs
        xrefs = xrefs = pull_xref(accession, from_database = 'P_REFSEQ_AC', to_database = 'ACC', xref_db_name = 'swissprot')
        return xrefs

def version_search(row):
    """
    Finds specific Uniprot version RefSeq cross-reference
    maps to if specific mapping exists
    """
    search = re.search('\[(.*?)\]', row['RefSeq'])
    if search:
        return search.group(1)
    return None

def expand_xref(xref):
    """
    Expands dataframe of Uniprot-RefSeq cross references
    """
    xref.dropna(subset = ['Cross-reference (RefSeq)'], inplace = True)
    xref['RefSeq']= xref['Cross-reference (RefSeq)'].str[:-1]
    
    xref_expanded = pd.DataFrame(xref['RefSeq'].str.split(';').tolist(), index=xref['Entry']).stack()
    xref_expanded = pd.DataFrame(xref_expanded).reset_index()
    
    xref_expanded.rename(columns = {0 : 'RefSeq'}, inplace = True)
    
    xref_expanded['Entry Full'] = xref_expanded.apply(version_search, axis = 1)
    xref_expanded['Entry Version'] = xref_expanded['Entry Full'].str.split('-').str[1]
    xref_expanded['RefSeq'] = xref_expanded['RefSeq'].apply(lambda x : re.sub("[\(\[].*?[\)\]]", "", x))
    xref_expanded['RefSeq ID'] = xref_expanded['RefSeq'].str.split('.').str[0]
    xref_expanded['RefSeq Verion'] = xref_expanded['RefSeq'].str.split('.').str[1]
    
    xref_expanded.drop(columns = ['level_1'], inplace = True)
    
    xref_expanded.rename(columns = {
        'Entry' : 'Uniprot ID',
        'Entry Full' : 'Uniprot ID Full',
        'Entry Version' : 'Uniprot ID Version',
        'RefSeq ID' : 'RefSeq ID',
        'RefSeq' : 'RefSeq Full',
        'RefSeq Version' : 'RefSeq Version',
        
    }, inplace = True)
    return xref_expanded

def full_xref_uniprot_refseq(reviewed = True):
    """
    All cross-references between Uniprot and RefSeq pulled from Uniprot database
    If a Uniprot Id does not map to a RefSeq id it is not included in the final dataframe
    
    Parameters
    ----------
    reviewed : bool
        True if only including reviewed Uniprot IDs
    
    Returns
    -------
    full_xref : pandas df
        Uniprot ID           Uniprot ID (P00533)
        RefSeq Full          Full RefSeq id (NP_005219.2)
        Uniprot ID Full      Full Uniprot ID with version (P00533-1). Only included if 
                                  specific RefSeq-Uniprot Full mapping
        Uniprot ID Version   Uniprot reference version
        RefSeq ID            RefSeq id without version (NP005219)
        RefSeq Version       RefSeq version 
    """
    url = 'https://www.uniprot.org/uniprot/'
    rev = 'yes' if reviewed else  'no'
        
    params = {
        'query' : f'reviewed:{rev} ',
        'columns' : 'id,database(RefSeq)',
        'format' : 'tab',
        'limit' : 1000
    }
    response = requests.get(url, params)
    
    if response.ok and response.text:
        df = pd.read_csv(StringIO(response.text), sep = '\t')
        full_xref = expand_xref(df)
        return full_xref
    return None

# import xml.dom.minidom as xml
# from xml.parsers.expat import ExpatError
# import urllib2, httplib
# from ptmscout.config import settings
# from ptmscout.utils.decorators import rate_limit
# import logging
# import socket

# log = logging.getLogger('ptmscout')

# class PICRParser(object):

#     def __init__(self, picrxml):
#         self.dom = xml.parseString(picrxml)
        
#         cross_refs = self.dom.getElementsByTagName('identicalCrossReferences')
#         self.references = []
        
#         for cross_ref in cross_refs:
#             self.parseCrossRef(cross_ref)
    
#     def __get_node_text(self, node):
#         txt = ""
#         for c in node.childNodes:
#             txt += str(c.nodeValue) 
#         return txt
    
#     def parseCrossRef(self, cross_ref):
#         dbName = self.__get_node_text(cross_ref.getElementsByTagName('databaseName')[0])
#         accession = self.__get_node_text(cross_ref.getElementsByTagName('accession')[0])
#         version = self.__get_node_text(cross_ref.getElementsByTagName('accessionVersion')[0])
        
#         self.references.append((dbName.lower(), accession))
#         self.references.append((dbName.lower(), "%s.%s" % (accession, version)))

# picr_accession_query_url = "http://www.ebi.ac.uk/Tools/picr/rest/getUPIForAccession?accession=%s&taxid=%d&database=REFSEQ&database=SWISSPROT"
# picr_sequence_query_url  = "http://www.ebi.ac.uk/Tools/picr/rest/getUPIForSequence?sequence=%s&taxid=%d&database=REFSEQ&database=SWISSPROT"
# PICR_QUERY_TIMEOUT=15
# PICR_QUERY_RETRIES=3


# class PICRError(Exception):
#     pass

# @rate_limit(rate=3)
# def get_picr(accession, taxon_id):
#     if settings.DISABLE_PICR:
#         return []

#     i = 0
#     while i  < PICR_QUERY_RETRIES:
#         i+=1
#         try:
#             result = urllib2.urlopen(picr_accession_query_url % (accession, taxon_id), timeout=PICR_QUERY_TIMEOUT)
#             xml_result = result.read()
#             if xml_result == '':
#                 return []

#             return PICRParser(xml_result).references
#         except ExpatError, e:
#             log.warning("ExpatError '%s' when querying PICR, retrying (%d / %d)", str(e), i, PICR_QUERY_RETRIES)
#         except urllib2.HTTPError, e:
#             log.warning("HTTPError '%s' when querying PICR, retrying (%d / %d)", str(e), i, PICR_QUERY_RETRIES)
#         except httplib.BadStatusLine:
#             log.warning("Got bad status line from PICR, retrying (%d / %d)", i, PICR_QUERY_RETRIES)
#         except socket.timeout, e:
#             log.warning("PICR Timeout reached, retrying (%d / %d)", i, PICR_QUERY_RETRIES)

#     raise PICRError()
