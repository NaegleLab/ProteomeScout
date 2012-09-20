from tests.DBTestCase import DBTestCase
import ptmscout.database.protein as protein

class ExpressionTestCase(DBTestCase):
    def test_gene_expression_should_ref_appropriate_data(self):
        prot = protein.getProteinById(35546)
        probe_ids = [ probe.probeset_id for probe in prot.expression_probes ]
        
        self.assertEqual(['203838_s_at', '203839_s_at', '216439_at'], probe_ids)
        
        probe = prot.expression_probes[0]
        
        self.assertEqual('203838_s_at', probe.probeset_id)
        self.assertEqual(162, len(probe.samples))
        
        human_samples = [ s for s in probe.samples if s.collection.name == 'Human' ]
        tounge_sample = [ s for s in probe.samples if s.tissue.name == 'Tongue' ]
        
        self.assertEqual(79, len(human_samples))
        self.assertEqual(1, len(tounge_sample))