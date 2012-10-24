import unittest
from ptmscout.utils.protein_utils import get_accession_type

class ProteinUtilsTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)

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