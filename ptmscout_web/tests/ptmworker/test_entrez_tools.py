from tests.PTMScoutTestCase import IntegrationTestCase
from ptmscout.config import settings
from Bio import Entrez
from ptmworker import entrez_tools, pfam_tools
import mock
from geeneus import Proteome

class EntrezQueryTestCase(IntegrationTestCase):

    def test_get_pubmed_record(self):
        record = entrez_tools.get_pubmed_record_by_id(12230038)
        
        self.assertEqual(settings.adminEmail, Entrez.email)
        
        self.assertEqual('The Bio* toolkits--a brief overview.', record['TI'])
        self.assertEqual('296-302', record['PG'])
        self.assertEqual('2002 Sep', record['DP'])

    def test_get_alignment_scores(self):
        seq1 = "MASWESRKLLLLLWKNFTLKRRKFGTLVSEIVLVLLLSIVLLTTRHLLSIKKIEALYFPDQPISTVPSFFR"
        seq2 = "KSDFMASWESRKLLLLLWKNFTLKRRKFGTLLTTRHLLSIKKIEALYFPDQPISTVPSFFRABDFGALYFPDQP"
        
        inserted, deleted = entrez_tools.get_alignment_scores(seq1, seq2)
        
        self.assertEqual(17, inserted)
        self.assertEqual(14, deleted)
        
    def test_map_domain_to_sequence_should_map_if_matched(self):
        seq1 = "MASWESRKLLLLLWKNFTLKRRKFGTLVSEIVLVLLLSIVLLTTRHLLSIKKIEALYFPDQPISTVPSFFR"
        seq2 = "MASWESRKLLLLLWKNFTLKRRKFGTLLTTRHLLSIKKIEALYFPDQPISTVPSFFR"
        domain_seq = "SIKKIEALYFPDQPI"
        
        pfd = pfam_tools.PFamDomain()
        pfd.accession = 'A034KA'
        pfd.class_ = 'Domain'
        pfd.significant = 1.3
        pfd.start = seq1.find(domain_seq)
        pfd.stop = pfd.start + len(domain_seq) - 1
        pfd.p_value = 0.05
        pfd.release = '1.2'
        pfd.label = 'Some domain'
        
        pfd2 = entrez_tools.map_domain_to_sequence(seq1, pfd, seq2)
        
        self.assertEqual(pfd.accession, pfd2.accession)
        self.assertEqual(pfd.class_, pfd2.class_)
        self.assertEqual(pfd.significant, pfd2.significant)
        self.assertEqual(pfd.p_value, pfd2.p_value)
        self.assertEqual(pfd.release, pfd2.release)
        self.assertEqual(pfd.label, pfd2.label)
        
        self.assertEqual(seq2.find(domain_seq), pfd2.start)
        self.assertEqual(seq2.find(domain_seq) + len(domain_seq) - 1, pfd2.stop)
        
    def test_map_domain_to_sequence_should_fail_if_not_matched(self):
        seq1 = "MASWESRKLLLLLWKNFTLKRRKFGTLVSEIVLVLLLSIVLLTTRHLLSIKKIEALYFPDQPISTVPSFFR"
        seq2 = "MASWESRKLLLLLWKNFTLKRRKFGTLLTTRHLLSIKKIEALYFPDQPISTVPSFFR"
        domain_seq = "GTLVSEIVLVLLLSIVLLT"
        
        pfd = pfam_tools.PFamDomain()
        pfd.accession = 'A034KA'
        pfd.class_ = 'Domain'
        pfd.significant = 1.3
        pfd.start = seq1.find(domain_seq)
        pfd.stop = pfd.start + len(domain_seq) - 1
        pfd.p_value = 0.05
        pfd.release = '1.2'
        pfd.label = 'Some domain'
        
        pfd2 = entrez_tools.map_domain_to_sequence(seq1, pfd, seq2)
        self.assertEqual(None, pfd2)
        

    def test_get_protein_information_should_raise_entrez_error(self):
        mock_pm = mock.Mock(spec=Proteome.ProteinManager)
        mock_pm.get_protein_sequence.return_value = ""

        try:
            entrez_tools.get_protein_information(mock_pm, 'P50914')
        except entrez_tools.EntrezError, e:
            self.assertEqual('P50914', e.acc)
        else:
            self.fail("Expected EntrezError")

    def test_get_isoform_map(self):
        accs, isoform_map = entrez_tools.get_isoform_map(['A6NC57', 'Q91ZU6-2', 'Q91ZU6-3', 'Q969I3-2'])

        self.assertEqual(set(['A6NC57', 'Q91ZU6', 'Q969I3']), set(accs))
        self.assertEqual({'Q969I3-2': 'Q969I3', 'Q91ZU6-2':'Q91ZU6', 'Q91ZU6-3': 'Q91ZU6'}, isoform_map)
        
    def test_get_protein_information(self):
        pm = Proteome.ProteinManager(settings.adminEmail)
        
        name, gene, taxonomy, species, _accessions, domains, seq, isoforms = entrez_tools.get_protein_information(pm, 'A6NC57')
        
        self.assertEqual('Ankyrin repeat domain-containing protein 62', name) 
        self.assertEqual(None, gene)
       
        taxons = [u'Eukaryota', u'Metazoa', u'Chordata', u'Craniata',
                u'Vertebrata', u'Euteleostomi', u'Mammalia', u'Eutheria',
                u'Euarchontoglires', u'Primates', u'Haplorrhini',
                u'Catarrhini', u'Hominidae', u'Homo']

        self.assertEqual('Homo sapiens', species)
        self.assertEqual([t.lower() for t in taxons], taxonomy)
        self.assertEqual('MEVRGSFLAACRRRMATWRKNRDKDGFSNPGYRVRQKDLGMIHKAAIAGDVNKVMESILLRLNDLNDRDKKNRTALLLACAHGRPGVVADLVARKCQLNLTDSENRTALIKAVQCQEEVCASILLEHGANPNVRDMYGNTALHYAIDNENISMARKLLAYGADIEARSQDGHTSLLLAVNRKKEQMVAFLLKKKPDLTAIDNFGRTALILAARNGSTSVVYQLLQHNIDVFCQDISGWTAEDYAVASKFQAIRGMISEYKANKRCKSLQNSNSEQDLEMTSEGEQERLEGCESSQPQVEEKMKKCRNKKMEVSRNVHADDSDNYNDDVDELIHKIKNRKPDNHQSPGKENGEFDRLARKTSNEKSKVKSQIYFTDDLNDISGSSEKTSEDDELPYSDDENFMLLIEQSGMECKDFVSLSKSKNATAACGRSIEDQKCYCERLKVKFQKMKNNISVLQKVLSETDKTKSQSEHQNLQGKKKLCNLRFILQQQEEERIKAEELYEKDIEELKIMEEQYRTQTEVKKQSKLTLKSLEVELKTVRSNSNQNFHTHERERDLWQENHLMRDEIARLRLEIDTIKHQNQETENKYFKDIEIIKENNEDLEKTLKRNEEALTKTITRYSKELNVLMDENTMLNSELQKEKQSMSRLETEMESYRCRLAAALCDHDQRQSSKRDLQLAFQSTVNEWCHLQEDTNSHIQILSQQLSKAESTSSGLETELHYEREALKEKTLHIEHMQGVLSRTQRRLEDIEHMYQNDQPILEKYVRKQQSVEDGLFQLQSQNLLYQQQCNDARKKADNQEKTIINIQVKCEDTVEKLQAECRKLEENNKGLMKECTLLKERQCQYEKEKEEREVVRRQLQREVDDALNKQLLLEAMLEISSERRINLEDEAQSLKKKLGQMRSQVCMKLSMSTVTL', seq) 
        
        parsed_domains = set([])
        for d in domains:
            parsed_domains.add( (d.label, d.start, d.stop) )
        
        self.assertIn(('Ank_2', 76, 167), parsed_domains)
        self.assertIn(('Ank_2', 142, 234), parsed_domains)
#        self.assertIn(('Ank_2', 109, 201), parsed_domains)
#        self.assertIn(('Ank_2', 93, 168), parsed_domains)
#        self.assertIn(('Ank_2', 42, 135), parsed_domains)
#        self.assertIn(('Ank_2', 175, 267), parsed_domains)
#        self.assertIn(('Pfam-B_116174', 445, 915), parsed_domains)
#        self.assertIn(('Pfam-B_3422', 342, 443), parsed_domains)

#        self.assertEqual('MEVRGSFLAACRRRMATWRKNRDKDGFSNPGYRVRQKDLGMIHKAAIAGDVNKVMESILLRLNDLNDRDKKNRTALLLACAHGRPGVVADLVARKCQLNLTDSENRTALIKAVQCQEEVCASILLEHGANPNVRDMYGNTALHYAIDNENISMARKLLAYGADIEARSQDGHTSLLLAVNRKKEQMVAFLLKKKPDLTAIDNFGRTALILAARNGSTSVVYQLLQHNIDVFCQDISGWTAEDYAVASKFQAIRGMISEYKANKRCKSLQNSNSEQDLEMTSEGEQERLEGCESSQPQVEEKMKKCRNKKMEVSRNVHADDSDNYNDDVDELIHKIKNRKPDNHQSPGKENGEFDR', isoforms['3'])
#        self.assertEqual('MTSEGEQERLEGCESSQPQVEEKMKKCRNKKMEVSRNVHADDSDNYNDDVDELIHKIKNRKPDNHQSPGKENGEFDRLARKTSNEKSKVKSQIYFTDDLNDISGSSEKTSEDDELPYSDDENFMLLIEQSGMECKDFVSLSKSKNATAACGRSIEDQKCYCERLKVKFQKMKNNISVLQKVLSETDKTKSQSEHQNLQGKKKLCNLRFILQQQEEERIKAEELYEKDIEELKIMEEQYRTQTEVKKQSKLTLKSLEVELKTVRSNSNQNFHTHERERDLWQENHLMRDEIARLRLEIDTIKHQNQETENKYFKDIEIIKENNEDLEKTLKRNEEALTKTITRYSKELNVLMDENTMLNSELQKEKQSMSRLETEMESYRCRLAAALCDHDQRQSSKRDLQLAFQSTVNEWCHLQEDTNSHIQILSQQLSKAESTSSGLETELHYEREALKEKTLHIEHMQGVLSRTQRRLEDIEHMYQNDQPILEKYVRKQQSVEDGLFQLQSQNLLYQQQCNDARKKADNQEKTIINIQVKCEDTVEKLQAECRKLEENNKGLMKECTLLKERQCQYEKEKEEREVVRRQLQREVDDALNKQLLLEAMLEISSERRINLEDEAQSLKKKLGQMRSQVCMKLSMSTVTL', isoforms['2'])
