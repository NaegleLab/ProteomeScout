# Skip
from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker.entrez_query import get_pubmed_record_by_id, get_protein_records_by_accession
from ptmscout.config import settings
from Bio import Entrez

class EntrezQueryTestCase(IntegrationTestCase):
        
    def test_get_pubmed_record(self):
        record = get_pubmed_record_by_id(12230038)
        
        self.assertEqual(settings.adminEmail, Entrez.email)
        
        self.assertEqual('The Bio* toolkits--a brief overview.', record['TI'])
        self.assertEqual('296-302', record['PG'])
        self.assertEqual('2002 Sep', record['DP'])

    def test_get_protein_records(self):
        accessions = ['XP_619639', 'P07197', 'Q14160', 'O15320', 'Q8N9T8', 'P51858', 'Q9H6F5']
        records = get_protein_records_by_accession(accessions)
        
        self.assertEqual(7, len(records))
        for acc in accessions:
            self.assertEqual(1, len(records[acc]))
