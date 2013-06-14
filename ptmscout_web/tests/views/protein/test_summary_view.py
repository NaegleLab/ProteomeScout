from ptmscout.config import strings
from tests.PTMScoutTestCase import IntegrationTestCase

class TestProteinSummaryView(IntegrationTestCase):
    def test_integration(self):
        result = self.ptmscoutapp.get('/proteins/35546/summary')

        result.mustcontain(strings.protein_summary_page_title)
        result.mustcontain("activated p21cdc42Hs kinase [Homo sapiens]")
        result.mustcontain("gi|8922075")
        result.mustcontain("NP_005772")
        result.mustcontain("ACK")
        result.mustcontain("http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;amp;id=gi|8922075")
        result.mustcontain("http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;amp;id=NP_005772")


