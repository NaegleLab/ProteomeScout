from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker.helpers import uniprot_tools
from ptmscout.config import settings
import os

class TestUniprotQuery(IntegrationTestCase):
#'annotations', 'dbxrefs', 'description', 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'upper'

    def test_uniprot_region_parse(self):
        result_xml = open(os.path.join(settings.ptmscout_path, 'tests','ptmworker','helpers','uniprot_result.xml'), 'r')
        rval = uniprot_tools.handle_result(result_xml)

        pr = rval['Q91ZU6']

        parsed_features = []
        for f in pr.features:
            parsed_features.append((f.type, f.name, f.source, f.start, f.stop))

        self.assertEqual(('domain', 'Actin-binding', 'uniprot', 35, 252,), parsed_features[0])
        self.assertEqual(('domain', 'SH3', 'uniprot', 889, 941,), parsed_features[4])
        self.assertEqual(('repeat', 'Plectin 5', 'uniprot', 1848, 1889,), parsed_features[9])

    def test_uniprot_handle_result(self):
        result_xml = open(os.path.join(settings.ptmscout_path, 'tests','ptmworker','helpers','uniprot_result.xml'), 'r')
        rval = uniprot_tools.handle_result(result_xml)

        pr = rval['Q91ZU6']
        self.assertEqual('DYST_MOUSE', pr.locus)
        self.assertEqual('Q91ZU6', pr.query_accession)
        self.assertEqual(12, len( pr.other_accessions ))
        self.assertEqual('Mus musculus', pr.species)
        self.assertEqual(0, len(pr.mutations))
        self.assertEqual(None, pr.host_organism)

    def test_uniprot_query_should_retain_strain_or_isolate(self):
        result = uniprot_tools.get_uniprot_records(['P75471'])
        pr = result['P75471']

        self.assertEqual( 'Cytadherence high molecular weight protein 2', pr.name )
        self.assertEqual( 'hmw2', pr.gene )
        self.assertEqual( 'MNDTDKKFPLQPVYDTGFDD', pr.sequence[:20] )
        self.assertEqual( 'Mycoplasma pneumoniae (strain ATCC 29342 / M129)', pr.species )

    def test_uniprot_get_mutants(self):
        result = uniprot_tools.get_uniprot_records(['A6NC57'])

        pr = result['A6NC57']

        self.assertEqual(4, len(pr.mutations))
        m_test_list = set()
        for m in pr.mutations:
            self.assertEqual('Substitution (single)', m.mutationType)
            self.assertEqual('A6NC57', m.acc_id)
            m_test_list.add(( m.location, m.original, m.mutant ))

        expected_muts = set([(188, 'A', 'S'),
                             (265, 'C', 'R'),
                             (406, 'E', 'K'),
                             (613, 'A', 'T')])

        self.assertEqual(expected_muts, m_test_list)

    def test_get_mutations_with_missing_sequences(self):
        results = uniprot_tools.get_uniprot_records(['P50914', 'Q6KC79'])

        self.assertEqual( 23, len(results['Q6KC79'].mutations) )
        self.assertEqual( 1, len(results['P50914'].mutations) )

    def test_uniprot_get_viral_protein(self):
        result = uniprot_tools.get_uniprot_records(['P03264','P03070'])
        pr = result['P03264']
        self.assertEqual('Homo sapiens', pr.host_organism)

        pr = result['P03070']
        self.assertEqual('Macaca', pr.host_organism)

    def test_uniprot_get_protein_by_accessions(self):
        accs = ['A0K45A', 'Q96P48-6', 'Q91ZU6-5', 'Q91ZU6-3', 'Q969I3-2', 'B1WC86', 'A1L112', 'B2RXC1', 'A6NF89', 'A6NH21']
        result = uniprot_tools.get_uniprot_records(accs)

        self.assertLess(set(accs[1:] + ['Q96P48']), set(result.keys()))
        self.assertEqual(23, len(result))

        pr = result['Q96P48']

        self.assertEqual(358, pr.mutations[0].location)
        self.assertEqual('Q', pr.mutations[0].mutant)
        self.assertEqual(1047, pr.mutations[1].location)
        self.assertEqual('E', pr.mutations[1].mutant)

