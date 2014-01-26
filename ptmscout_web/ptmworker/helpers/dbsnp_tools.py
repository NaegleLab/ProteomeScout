from Bio import Entrez
from ptmscout.config import settings
from ptmworker.helpers import upload_helpers
from xml.dom import minidom, Node
import re

def getText(pNode):
    if pNode.nodeType == Node.TEXT_NODE:
        return pNode.nodeValue

    ts = ""
    for cNode in pNode.childNodes:
        ts += getText(cNode)
    return ts

def getChildNodesByTag(tag, entry):
    cNodes = []
    for cNode in entry.childNodes:
        if cNode.localName.lower() == tag.lower():
            cNodes.append(cNode)

    return cNodes

def parseRSEntry(entry):
    rsId = "rs%s" % (entry.attributes['rsId'].nodeValue)

    merged_ids = []
    mhNodes = getChildNodesByTag('MergeHistory', entry)
    for mhNode in mhNodes:
        merged_rsId = "rs%s" % (mhNode.attributes['rsId'].nodeValue)
        merged_ids.append(merged_rsId)

    clinicalPhenotypes = []

    phenotypes = getChildNodesByTag('Phenotype', entry)
    for pNode in phenotypes:
        csNodes = getChildNodesByTag('ClinicalSignificance', pNode)
        for csNode in csNodes:
            clinicalPhenotypes.append(getText(csNode))

    return rsId, merged_ids, ";".join(clinicalPhenotypes)

def parseRecords(handle):
    result = {}
    res = minidom.parseString(handle.read())

    rsEntries = res.getElementsByTagName('Rs')
    for rsEntry in rsEntries:
        rsId, mergedIds, clinicalPhenotype = parseRSEntry(rsEntry)

        result[rsId] = clinicalPhenotype
        for merged_rsId in mergedIds:
            result[merged_rsId] = clinicalPhenotype
    return result

def get_variant_classification(dbsnp_ids):
    Entrez.email = settings.adminEmail
    handle = Entrez.efetch(db="snp", id=dbsnp_ids, retmode="xml")

    return parseRecords(handle)
