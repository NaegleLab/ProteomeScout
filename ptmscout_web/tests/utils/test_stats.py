import unittest
from ptmscout.utils import stats

class TestDatasetExplorerView(unittest.TestCase):
    
    def test_prime_sieve(self):
        self.assertEqual(168, len( stats.computeSeive(1000) ))
    
    def test_fisher_and_significance(self):
        evalue = stats.fisher(2,7,8,4)
        self.assertAlmostEqual(0.0505, evalue, 4)
        pvalue = stats.fisher_pvalue(evalue, 2,7,8,4)
        self.assertAlmostEqual(0.0805, pvalue, 4)
        
    def test_fisher_large_numbers(self):
        stats.fisher_init(50000)
        
        evalue = stats.fisher(12, 200, 400, 30000)
        self.assertAlmostEqual(2.692e-5, evalue, 8)
        
        pvalue = stats.fisher_pvalue(evalue, 12, 200, 400, 30000)
        self.assertAlmostEqual(3.371e-5, pvalue, 8)