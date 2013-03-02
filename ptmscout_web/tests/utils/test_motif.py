from tests.PTMScoutTestCase import IntegrationTestCase
from ptmscout.utils import motif
from ptmscout.database import modifications

class TestMotifWithTestDB(IntegrationTestCase):
    def test_calculate_motif_enrichment_should_calculate_enrichment(self):
        foreground_msids = set([399217,399218,399219,399229,399230,399228])
        background = modifications.getMeasuredPeptidesByExperiment(28)
        foreground = [ms for ms in background if ms.id in foreground_msids]

        test_cnt, results = motif.calculate_motif_enrichment(foreground, background)
        
        self.assertEqual(38, test_cnt)
        
        exp_results = [
                ('......EyO......', (5, 7), (6, 73), 8.18e-06), 
                ('......-y......+', (5, 7), (6, 73), 8.18e-06),   
                ('....O.EyO.....+', (4, 7), (4, 73), 3.22e-05), 
                ('.T....EyO......', (4, 7), (5, 73), 0.000155), 
                ('......-y.......', (6, 7), (15, 73), 0.000182), 
                ('HTGFLTEyVATRWYR', (3, 7), (3, 73), 0.000563), 
                ('..L..P-y.....PK', (2, 7), (2, 73), 0.00799), 
            ]
        self.assertEqual(exp_results, results)