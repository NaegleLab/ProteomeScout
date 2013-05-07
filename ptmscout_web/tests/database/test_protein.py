from tests.DBTestCase import DBTestCase
from ptmscout.database import protein
from ptmscout.database.protein import NoSuchProtein, GeneOntology
from sqlalchemy.exc import IntegrityError

class ProteinTest(DBTestCase):

    def test_hasPrediction_should_verify_prediction(self):
        p = protein.Protein()
        
        s1 = protein.ProteinScansite()
        s1.source = 'scansite_bind'
        s1.value = 'SRC'
        s1.site_pos = 800
        
        s2 = protein.ProteinScansite()
        s2.source = 'scansite_kinase'
        s2.value = 'Pkin_S'
        s2.site_pos = 800
        
        s3 = protein.ProteinScansite()
        s3.source = 'scansite'
        s3.value = 'Something'
        s3.site_pos = 900
        
        p.scansite = [s1, s2, s3]
        
        self.assertTrue( p.hasPrediction('scansite_bind', 'SRC', 800) )
        self.assertFalse( p.hasPrediction('scansite_bind', 'SRC', 801) )
        self.assertFalse( p.hasPrediction('scansite_bind', 'SRC1', 800) )
        self.assertFalse( p.hasPrediction('scansite_kinase', 'SRC', 800) )

    def test_gene_ontology_hasChild_should_find_child_terms(self):
        go = GeneOntology()
        go.aspect = 'F'
        go.GO = 'GO:00000001'
        go.term = "SomeTermOfGO"
        go.version = "1.2"
        
        c1 = GeneOntology()
        c1.aspect = 'F'
        c1.GO = 'GO:00000002'
        c1.term = "SomeTermOfGO"
        c1.version = "1.2"

        c2 = GeneOntology()
        c2.aspect = 'F'
        c2.GO = 'GO:00000003'
        c2.term = "SomeTermOfGO"
        c2.version = "1.2"

        go.children = [c1,c2]

        self.assertTrue(go.hasChild('GO:00000003'))
        self.assertFalse(go.hasChild('GO:00000001'))
        self.assertFalse(go.hasChild('GO:00000004'))

    def test_protein_should_be_gettable(self):
        prot = protein.getProteinById(14496)
        
        self.assertEqual(14496, prot.id)
        self.assertEqual('MAAVILESIFLKRSQQKKKTSPLNFKKRLFLLTVHKLSYYEYDFERGRRGSKKGSIDVEKITCVETVVPEKNPPPERQIPRRGEESSEMEQISIIERFPYPFQVVYDEGPLYVFSPTEELRKRWIHQLKNVIRYNSDLVQKYHPCFWIDGQYLCCSQTAKNAMGCQILENRNGSLKPGSSHRKTKKPLPPTPEEDQILKKPLPPEPAAAPVSTSELKKVVALYDYMPMNANDLQLRKGDEYFILEESNLPWWRARDKNGQEGYIPSNYVTEAEDSIEMYEWYSKHMTRSQAEQLLKQEGKEGGFIVRDSSKAGKYTVSVFAKSTGDPQGVIRHYVVCSTPQSQYYLAEKHLFSTIPELINYHQHNSAGLISRLKYPVSQQNKNAPSTAGLGYGSWEIDPKDLTFLKELGTGQFGVVKYGKWRGQYDVAIKMIKEGSMSEDEFIEEAKVMMNLSHEKLVQLYGVCTKQRPIFIITEYMANGCLLNYLREMRHRFQTQQLLEMCKDVCEAMEYLESKQFLHRDLAARNCLVNDQGVVKVSDFGLSRYVLDDEYTSSVGSKFPVRWSPPEVLMYSKFSSKSDIWAFGVLMWEIYSLGKMPYERFTNSETAEHIAQGLRLYRPHLASEKVYTIMYSCWHEKADERPTFKILLSNILDVMDEES', prot.sequence)
        self.assertEqual('homo sapiens', prot.species.name)
        self.assertEqual('BTK', prot.acc_gene)
        self.assertEqual('Tyrosine-protein kinase BTK; AltName: Full=Bruton tyrosine kinase; AltName: Full=Agamm', prot.name)
         
    def test_non_existant_protein_should_raise_no_such_protein(self):
        try:
            pid = 10000000
            protein.getProteinById(pid)
        except NoSuchProtein, nsp:
            self.assertEqual(pid, nsp.prot)
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception: NoSuchProtein")
        
    def test_proteins_should_be_associated_with_correct_accession(self):
        prot = protein.getProteinById(1)
        
        self.assertEqual('USP24', prot.acc_gene)
        
        acc_ids = set([ acc.value for acc in prot.accessions ])
        self.assertEqual(
                set(['IPI00398505','NP_056121','Q9UPU5','IPI00902614', 'KIAA1057']), 
                acc_ids)
        
    def test_proteins_should_be_associated_with_correct_GO_terms(self):
        prot = protein.getProteinById(100)
        
        self.assertEqual('Mpv17l2', prot.acc_gene)
        
        acc_ids = set([ entry.GO_term.term for entry in prot.GO_terms ])
        self.assertEqual(
                set(['mitochondrion','molecular_function','biological_process']), 
                acc_ids)

    def test_GO_terms_must_be_unique(self):
        prot = protein.getProteinById(1)
        
        go_term1 = GeneOntology()
        go_term1.aspect = 'F'
        go_term1.GO = 'GO:unique'
        go_term1.term = "SomeTermOfGO"
        go_term1.version = "1.2"
        
        go_term2 = GeneOntology()
        go_term2.aspect = 'F'
        go_term2.GO = 'GO:unique'
        go_term2.term = "SomeTermOfGO"
        go_term2.version = "1.2"
        
        prot.addGoTerm(go_term1)
        prot.addGoTerm(go_term2)
        
        try:
            prot.saveProtein()
        except IntegrityError:
            pass
        except Exception, e:
            self.fail("Unexpected exception occurred: " + str(e))
        else:
            self.fail("Expected IntegrityError on last DB operation")
   
    def test_searchProteins(self):
        cnt, prots = protein.searchProteins("ACK1")

        expected_ids = [10367, 17551, 19946, 35546]
        self.assertEqual(expected_ids, sorted( [p.id for p in prots] ))
        self.assertEqual(4, cnt)

    def test_searchProteins_filter_species(self):
        expected_ids = [35546, 10367]
        
        cnt, prots = protein.searchProteins("ACK1", "homo sapiens")
        
        self.assertEqual(sorted(expected_ids), sorted([ prot.id for prot in prots ]))
        self.assertEqual(2, cnt)

    def test_searchProteins_by_acc_gene(self):
        expected_ids = [1212, 13807]

        cnt, prots = protein.searchProteins("CLEC16A")

        self.assertEqual(2, cnt)
        self.assertEqual(expected_ids, sorted( [p.id for p in prots] ))

    def test_searchProteins_by_name(self):
        cnt, prots = protein.searchProteins("ubiquitin-protein ligase", "homo sapiens", includeNames=True)
        self.assertEqual(88, cnt)
        self.assertEqual(88, len(prots))

    def test_searchProteins_by_sequence(self):
        exp_ids = [7, 1212]
        cnt, prots = protein.searchProteins(sequence="HGKTSRNIHSLDHLKYLYHVLT")

        self.assertEqual(2, cnt)
        self.assertEqual(exp_ids, sorted( [ p.id for p in prots ] ))
 
    def test_getProteinsByGene(self):
        expected_ids = [10367, 17551, 19946, 35546]
        
        prots = protein.getProteinsByGene("ACK1")
        
        self.assertEqual(sorted(expected_ids), sorted([ prot.id for prot in prots ]))
    
    def test_getProteinsByGene_filter_species(self):
        expected_ids = [35546, 10367]
        
        prots = protein.getProteinsByGene("ACK1", "homo sapiens")
        
        self.assertEqual(sorted(expected_ids), sorted([ prot.id for prot in prots ]))
                
             
    
    def test_proteins_should_be_associated_with_correct_domains(self):
        prot = protein.getProteinById(100)
        
        self.assertEqual('Mpv17l2', prot.acc_gene)
        
        domains = set([ d.label for d in prot.domains ])
        self.assertEqual( set(['Mpv17_PMP22']), domains )
        
    def test_accessions_should_be_sorted_by_types(self):
        prot = protein.getProteinById(10367)
        
        types = [acc.type for acc in prot.accessions]
        
        self.assertEqual(sorted(types), types)
