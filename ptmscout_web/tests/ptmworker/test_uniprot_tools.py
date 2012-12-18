from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker import uniprot_tools

class TestUniprotQuery(IntegrationTestCase):
#'annotations', 'dbxrefs', 'description', 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'upper'

    def test_uniprot_handle_result(self):
        result_xml = open('uniprot_result.xml', 'r')
        uniprot_tools.handle_result(result_xml)

    def test_uniprot_get_protein_by_accessions(self):
        accs = ['A0K45A', 'Q91ZU6-5', 'Q91ZU6-3', 'Q969I3-2', 'B1WC86', 'A1L112', 'B2RXC1', 'A6NF89', 'A6NH21']
        result = uniprot_tools.get_uniprot_records(accs)

        self.assertTrue(set(accs[1:]) < set(result.keys()))
        self.assertEqual(14, len(result))
