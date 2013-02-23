import unittest
from ptmscout.utils.protein_utils import get_accession_type, create_sequence_profile
from tests.views.mocking import createMockProtein, createMockMeasurement,\
    createMockPeptide, createMockPTM, createMockPeptideModification
import math

class ProteinUtilsTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def test_get_sequence_profile_with_empty_set(self):
        result = create_sequence_profile([])
        
        self.assertEqual({'total':0, 'frequencies':[]}, result)
    
    def test_get_sequence_profile(self):
        p1 = createMockProtein()
        
        m1 = createMockMeasurement(p1.id, 1)
        m2 = createMockMeasurement(p1.id, 1)
        
        pep1 = createMockPeptide(p1.id)
        pep2 = createMockPeptide(p1.id)
        pep3 = createMockPeptide(p1.id)
        
        mod = createMockPTM()
        
        createMockPeptideModification(m1, pep1, mod)
        createMockPeptideModification(m1, pep2, mod)
        createMockPeptideModification(m2, pep3, mod)
        
        pep1.pep_aligned='LKKVVALyDYMPMNA'
        pep2.pep_aligned=' SHWQQQsYLDSGIH'
        pep3.pep_aligned=' ATWTAQsLLGSGIP'
        
        mods = [m1,m2]
        result = create_sequence_profile(mods)
        
        en = 19 / (2 * math.log(2) * 3)
        R21 = math.log(20, 2) + ( (2/3.0) * math.log((2/3.0), 2) + (1/3.0) * math.log((1/3.0), 2) ) - en
        R111 = math.log(20, 2) + ( (1/3.0) * math.log((1/3.0), 2) + (1/3.0) * math.log((1/3.0), 2) + (1/3.0) * math.log((1/3.0), 2) ) - en
        
        expected_peps = []
        expected_peps.append({'R':R21, 'f':[('-', 2, 0),('L', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('A', 1, 0),('K', 1, 1),('S', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('H', 1, 0),('K', 1, 1),('T', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('W', 2, 0),('V', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('Q', 1, 0),('T', 1, 1),('V', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('A', 2, 0),('Q', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('Q', 2, 0),('L', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('S', 2, 0),('Y', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('Y', 1, 0),('D', 1, 1),('L',1, 2)]})
        expected_peps.append({'R':R21, 'f':[('L', 2, 0),('Y', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('M', 1, 0),('D', 1, 1),('G',1, 2)]})
        expected_peps.append({'R':R21, 'f':[('S', 2, 0),('P', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('G', 2, 0),('M', 1, 2)]})
        expected_peps.append({'R':R21, 'f':[('I', 2, 0),('N', 1, 2)]})
        expected_peps.append({'R':R111, 'f':[('A', 1, 0),('H', 1, 1),('P',1, 2)]})
        
        expected = {'total': 3, 'frequencies': expected_peps}
        
        self.assertEqual(expected['total'], result['total'])
        self.assertEqual( len(expected['frequencies']), len(result['frequencies']) )
        
        for i in xrange(0, 15):
            exp = expected['frequencies'][i]
            res = result['frequencies'][i]
            
            self.assertAlmostEqual(exp['R'], res['R'], 6)
            self.assertEqual(exp['f'], res['f'])    

    def test_get_accession_type(self):
        self.assertEqual('gi',          get_accession_type('gi|13245'))
        self.assertEqual('refseq',      get_accession_type('NP_001101107'))
        self.assertEqual('swissprot',   get_accession_type('Q9NQG7'))
        self.assertEqual('swissprot',   get_accession_type('NEBU_HUMAN'))
        self.assertEqual('genbank',     get_accession_type('CAD28549'))
        self.assertEqual('ipi',         get_accession_type('IPI00902614'))
        self.assertEqual('ensembl',    get_accession_type('ENSVPAP00000000821'))
        
if __name__ == '__main__':
    unittest.main()