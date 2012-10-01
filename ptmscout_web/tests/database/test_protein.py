from tests.DBTestCase import DBTestCase
from ptmscout.database import protein
from ptmscout.database.protein import NoSuchProtein, GeneOntology
from sqlalchemy.exc import IntegrityError

class ProteinTest(DBTestCase):
    def test_protein_should_be_gettable(self):
        prot = protein.getProteinById(14496)
        
        self.assertEqual(14496, prot.id)
        self.assertEqual('MAAVILESIFLKRSQQKKKTSPLNFKKRLFLLTVHKLSYYEYDFERGRRGSKKGSIDVEKITCVETVVPEKNPPPERQIPRRGEESSEMEQISIIERFPYPFQVVYDEGPLYVFSPTEELRKRWIHQLKNVIRYNSDLVQKYHPCFWIDGQYLCCSQTAKNAMGCQILENRNGSLKPGSSHRKTKKPLPPTPEEDQILKKPLPPEPAAAPVSTSELKKVVALYDYMPMNANDLQLRKGDEYFILEESNLPWWRARDKNGQEGYIPSNYVTEAEDSIEMYEWYSKHMTRSQAEQLLKQEGKEGGFIVRDSSKAGKYTVSVFAKSTGDPQGVIRHYVVCSTPQSQYYLAEKHLFSTIPELINYHQHNSAGLISRLKYPVSQQNKNAPSTAGLGYGSWEIDPKDLTFLKELGTGQFGVVKYGKWRGQYDVAIKMIKEGSMSEDEFIEEAKVMMNLSHEKLVQLYGVCTKQRPIFIITEYMANGCLLNYLREMRHRFQTQQLLEMCKDVCEAMEYLESKQFLHRDLAARNCLVNDQGVVKVSDFGLSRYVLDDEYTSSVGSKFPVRWSPPEVLMYSKFSSKSDIWAFGVLMWEIYSLGKMPYERFTNSETAEHIAQGLRLYRPHLASEKVYTIMYSCWHEKADERPTFKILLSNILDVMDEES', prot.sequence)
        self.assertEqual('homo sapiens', prot.species.name)
        self.assertEqual('BTK', prot.acc_gene)
        self.assertEqual('Tyrosine-protein kinase BTK; AltName: Full=Bruton tyrosine kinase; AltName: Full=Agamm', prot.name)
        self.assertEqual('12-2009', prot.date)
         
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
        
        acc_ids = set([ term.term for term in prot.GO_terms ])
        self.assertEqual(
                set(['mitochondrion','molecular_function','biological_process']), 
                acc_ids)

    def test_GO_terms_must_be_unique(self):
        prot = protein.getProteinById(1)
        
        go_term1 = GeneOntology()
        go_term1.aspect = 'F'
        go_term1.GO = 'GO:unique'
        go_term1.term = "SomeTermOfGO"
        go_term2 = GeneOntology()
        go_term2.aspect = 'F'
        go_term2.GO = 'GO:unique'
        go_term2.term = "SomeTermOfGO"
        
        prot.GO_terms.append(go_term1)
        prot.GO_terms.append(go_term2)
        prot.saveProtein()
        
        try:
            self.session.flush()
        except IntegrityError:
            pass
        except Exception, e:
            self.fail("Unexpected exception occurred: " + str(e))
        else:
            self.fail("Expected IntegrityError on last DB operation")
    
    def test_getProteinsByAccession(self):
        expected_ids = [27866, 35546, 10367, 30840, 19946, 17551, 35433]
        
        prots = protein.getProteinsByAccession(["ACK1"])
        
        self.assertEqual(sorted(expected_ids), sorted([ prot.id for prot in prots ]))
    
    def test_getProteinsByAccession_filter_species(self):
        expected_ids = [35546, 10367]
        
        prots = protein.getProteinsByAccession(["ACK1"], "homo sapiens")
        
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