from tests.DBTestCase import DBTestCase
from ptmscout.database import protein, mutations

class TestMutations(DBTestCase):
    def test_protein_has_mutation(self):
        p = protein.getProteinById(35546)
        m1 = mutations.Mutation('single', 134, 'A', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)

        m2 = mutations.Mutation('single', 156, 'G', 'D', 'ACK1_HUMAN',
                 'some annotations', p.id)


        test_m1 = mutations.Mutation('single', 190, 'L', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)
        test_m2 = mutations.Mutation('single', 134, 'L', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)
        test_m3 = mutations.Mutation('single', 156, 'L', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)

        self.assertFalse(p.hasMutation(test_m1))
        self.assertFalse(p.hasMutation(test_m2))
        self.assertFalse(p.hasMutation(test_m3))

        p.mutations.extend([m1, m2])

        self.assertFalse(p.hasMutation(test_m1))
        self.assertTrue(p.hasMutation(test_m2))
        self.assertFalse(p.hasMutation(test_m3))



    def test_equals_works(self):
        p = protein.getProteinById(35546)
        m1 = mutations.Mutation('single', 134, 'A', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)

        m2 = mutations.Mutation('single', 156, 'G', 'D', 'ACK1_HUMAN',
                 'some annotations', p.id)


        test_m1 = mutations.Mutation('single', 190, 'L', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)
        test_m2 = mutations.Mutation('single', 134, 'L', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)
        test_m3 = mutations.Mutation('single', 156, 'L', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)

        self.assertTrue( m1.equals(test_m2) )
        self.assertTrue( test_m2.equals(m1) )

        self.assertFalse( m1.equals(test_m1) )
        self.assertFalse( test_m1.equals(m1) )

        self.assertFalse( m2.equals(test_m2) )
        self.assertFalse( test_m2.equals(m2) )

        self.assertFalse( m1.equals(test_m3) )
        self.assertFalse( test_m3.equals(m1) )

        self.assertFalse( m2.equals(test_m3) )
        self.assertFalse( test_m3.equals(m2) )

    def test_is_consistent(self):
        p = protein.getProteinById(35546)

        # MQPEEGTGWL LELLSEVQLQ QYFLRLRDDL NVTRLSHFEY VKNEDLEKIG M...
        test_m1 = mutations.Mutation('single', 10, 'L', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)
        test_m2 = mutations.Mutation('single', 10, 'D', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)
        test_m3 = mutations.Mutation('single', 11, 'LELL', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)
        test_m4 = mutations.Mutation('single', 1, 'MQPQ', 'K', 'ACK1_HUMAN',
                 'some annotations', p.id)

        self.assertTrue( test_m1.consistent( p.sequence ) )
        self.assertFalse( test_m2.consistent( p.sequence ) )
        self.assertTrue( test_m3.consistent( p.sequence ) )
        self.assertFalse( test_m4.consistent( p.sequence ) )
