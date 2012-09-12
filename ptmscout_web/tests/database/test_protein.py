from tests.DBTestCase import DBTestCase
from ptmscout.database import protein
from ptmscout.database.protein import NoSuchProtein

class ProteinTest(DBTestCase):
    def test_protein_should_be_gettable(self):
        prot = protein.getProteinById(14496)
        
        self.assertEqual(14496, prot.id)
        self.assertEqual('MAAVILESIFLKRSQQKKKTSPLNFKKRLFLLTVHKLSYYEYDFERGRRGSKKGSIDVEKITCVETVVPEKNPPPERQIPRRGEESSEMEQISIIERFPYPFQVVYDEGPLYVFSPTEELRKRWIHQLKNVIRYNSDLVQKYHPCFWIDGQYLCCSQTAKNAMGCQILENRNGSLKPGSSHRKTKKPLPPTPEEDQILKKPLPPEPAAAPVSTSELKKVVALYDYMPMNANDLQLRKGDEYFILEESNLPWWRARDKNGQEGYIPSNYVTEAEDSIEMYEWYSKHMTRSQAEQLLKQEGKEGGFIVRDSSKAGKYTVSVFAKSTGDPQGVIRHYVVCSTPQSQYYLAEKHLFSTIPELINYHQHNSAGLISRLKYPVSQQNKNAPSTAGLGYGSWEIDPKDLTFLKELGTGQFGVVKYGKWRGQYDVAIKMIKEGSMSEDEFIEEAKVMMNLSHEKLVQLYGVCTKQRPIFIITEYMANGCLLNYLREMRHRFQTQQLLEMCKDVCEAMEYLESKQFLHRDLAARNCLVNDQGVVKVSDFGLSRYVLDDEYTSSVGSKFPVRWSPPEVLMYSKFSSKSDIWAFGVLMWEIYSLGKMPYERFTNSETAEHIAQGLRLYRPHLASEKVYTIMYSCWHEKADERPTFKILLSNILDVMDEES', prot.sequence)
        self.assertEqual('homo sapiens', prot.species)
        self.assertEqual('BTK', prot.acc_gene)
        self.assertEqual('Tyrosine-protein kinase BTK; AltName: Full=Bruton tyrosine kinase; AltName: Full=Agamm', prot.name)
        self.assertEqual('12-2009', prot.date)
         
    def test_non_existant_protein_should_raise_no_such_protein(self):
        try:
            pid = 10000000
            protein.getProteinById(pid)
        except NoSuchProtein, nsp:
            self.assertEqual(pid, nsp.pid)
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception: NoSuchProtein")
        
        