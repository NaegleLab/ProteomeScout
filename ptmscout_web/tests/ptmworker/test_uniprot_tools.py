from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker import uniprot_tools

class TestUniprotQuery(IntegrationTestCase):
#'annotations', 'dbxrefs', 'description', 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'upper'
    
    def test_uniprot_get_protein_isoforms(self):
        accs, isoform_map = uniprot_tools.get_isoform_map(['Q91ZU6-2', 'Q91ZU6-3', 'Q969I3-2'])

        uniprot_tools.get_protein_isoforms(isoform_map)
    
    def test_uniprot_get_protein_by_accessions(self):
        accs = ['Q91ZU6-1', 'Q91ZU6-3', 'Q969I3-2', 'A1L112', 'B1WC86', 'A1L112', 'B2RXC1', 'A6NF89', 'A6NH21']
        uniprot_tools.get_uniprot_records(accs)
        