from tests.PTMScoutTestCase import IntegrationTestCase
from ptmscout.config import settings
from ptmworker.helpers import dbsnp_tools
import mock

class DBSNPQueryTestCase(IntegrationTestCase):

    def test_get_dbsnp_id_list(self):
        snp_list = ["rs17336639", "rs17289589", "rs2227983"]

        result = dbsnp_tools.get_variant_classification(snp_list)

        self.assertEqual('', result['rs17336639'])
        self.assertEqual('', result['rs17289589'])
        self.assertEqual('', result['rs2227983'])

    def test_get_dbsnp_with_clinical_significance(self):
        snp_list= ["rs268", "rs328", "rs334"]

        result = dbsnp_tools.get_variant_classification(snp_list)

        self.assertEqual(3, len(result))
        self.assertEqual('pathogenic', result['rs268'])
        self.assertEqual('non-pathogenic', result['rs328'])
        self.assertEqual('pathogenic', result['rs334'])


    def test_get_merged_rsid(self):
        snp_list = ["rs17862494"]

        result = dbsnp_tools.get_variant_classification(snp_list)

        self.assertEqual(2, len(result))
        self.assertEqual('', result["rs17862494"])
