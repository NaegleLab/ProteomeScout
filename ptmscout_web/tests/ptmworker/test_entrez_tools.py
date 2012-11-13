from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker.entrez_query import get_pubmed_record_by_id
from ptmscout.config import settings
from Bio import Entrez

class EntrezQueryTestCase(IntegrationTestCase):
        
    def test_get_pubmed_record(self):
        record = get_pubmed_record_by_id(12230038)
        
        self.assertEqual(settings.adminEmail, Entrez.email)
        
        self.assertEqual('The Bio* toolkits--a brief overview.', record['TI'])
        self.assertEqual('296-302', record['PG'])
        self.assertEqual('2002 Sep', record['DP'])
