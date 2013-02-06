from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker.helpers import uniprot_tools

class TestUniprotQuery(IntegrationTestCase):
#'annotations', 'dbxrefs', 'description', 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'upper'

    def test_uniprot_handle_result(self):
        result_xml = open('uniprot_result.xml', 'r')
        uniprot_tools.handle_result(result_xml)

    def test_uniprot_query_should_retain_strain_or_isolate(self):
        result = uniprot_tools.get_uniprot_records(['P75471'])
        name, gene, taxons, species, host_organism, other_accessions, domains, mutations, seq = result['P75471']

        self.assertEqual( 'Cytadherence high molecular weight protein 2', name )
        self.assertEqual( 'hmw2', gene )
        self.assertEqual( 'MNDTDKKFPLQPVYDTGFDD', seq[:20] )
        self.assertEqual( 'Mycoplasma pneumoniae (strain ATCC 29342 / M129)', species )

    def test_uniprot_get_mutants(self):
        result = uniprot_tools.get_uniprot_records(['A6NC57'])

        name, gene, taxons, species, host_organism, other_accessions, domains, mutations, seq = result['A6NC57']

        self.assertEqual(4, len(mutations))
        m_test_list = set()
        for m in mutations:
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

        self.assertEqual( 23, len(results['Q6KC79'][7]) )
        self.assertEqual( 1, len(results['P50914'][7]) )

    def test_uniprot_get_viral_protein(self):
        result = uniprot_tools.get_uniprot_records(['P03264','P03070'])
        name, gene, taxons, species, host_organism, other_accessions, domains, mutations, seq = result['P03264']
        self.assertEqual('Homo sapiens', host_organism)

        name, gene, taxons, species, host_organism, other_accessions, domains, mutations, seq = result['P03070']
        self.assertEqual('Macaca', host_organism)

    def test_uniprot_get_protein_by_accessions(self):
        accs = ['A0K45A', 'Q96P48-6', 'Q91ZU6-5', 'Q91ZU6-3', 'Q969I3-2', 'B1WC86', 'A1L112', 'B2RXC1', 'A6NF89', 'A6NH21']
        result = uniprot_tools.get_uniprot_records(accs)

        self.assertLess(set(accs[1:] + ['Q96P48']), set(result.keys()))
        self.assertEqual(23, len(result))

        _n, _g, _t, _s, _h, _o, _d, mutations, _s = result['Q96P48']

        self.assertEqual(358, mutations[0].location)
        self.assertEqual('Q', mutations[0].mutant)
        self.assertEqual(1047, mutations[1].location)
        self.assertEqual('E', mutations[1].mutant)

