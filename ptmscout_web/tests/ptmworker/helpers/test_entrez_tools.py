from tests.PTMScoutTestCase import IntegrationTestCase
from ptmscout.config import settings
from Bio import Entrez
from ptmworker.helpers import entrez_tools, pfam_tools
import mock
from geeneus import Proteome

class EntrezQueryTestCase(IntegrationTestCase):

    def test_get_taxonomy_node(self):
        taxon = "Xanthomonas campestris pv. oryzae"

        taxonomy = entrez_tools.get_taxonomic_lineage(taxon)

        exp_lin = \
            [('cellular organisms', 131567), ('Bacteria', 2), ('Proteobacteria', 1224),
             ('Gammaproteobacteria', 1236), ('Xanthomonadales', 135614), ('Xanthomonadaceae', 32033),
             ('Xanthomonas', 338), ('Xanthomonas campestris', 339), ('Xanthomonas campestris pv. oryzae', 314227)]

        self.assertEqual(exp_lin, taxonomy)

    def test_get_pubmed_record_22373819(self):
        record = entrez_tools.get_pubmed_record_by_id(22373819)

        self.assertEqual("571-571", record['PG'])

    def test_get_pubmed_record_12230038(self):
        record = entrez_tools.get_pubmed_record_by_id(12230038)
        
        self.assertEqual(settings.adminEmail, Entrez.email)
        
        self.assertEqual('The Bio* toolkits--a brief overview.', record['TI'])
        self.assertEqual('296-302', record['PG'])
        self.assertEqual('2002 Sep', record['DP'])

    def test_get_viral_protein(self):
        result = entrez_tools.get_proteins_from_ncbi(['118734'])

        pr = result[0]['118734']

        self.assertEqual('DNA-binding protein; Short=DBP; AltName: Full=Early 2A protein', pr.name)


    def test_get_record_with_pipe_format_accessions(self):
        pm = Proteome.ProteinManager(settings.adminEmail, uniprotShortcut=False)
        pr = entrez_tools.get_protein_information(pm, 'CAI16470.1')
        self.assertEqual(3, len(pr.other_accessions))
        self.assertEqual(set([('gi', 'gi|55958964'), ('embl', 'CAI16470'), ('embl', 'CAI16470.1')]), set(pr.other_accessions))

        pr = entrez_tools.get_protein_information(pm, 'O75643')
        self.assertEqual(4, len(pr.other_accessions))
        self.assertEqual(set([('swissprot', 'O75643.2'), ('swissprot', 'U520_HUMAN'), ('gi', 'gi|56405304'), ('swissprot', 'O75643')]), set(pr.other_accessions))

        pr = entrez_tools.get_protein_information(pm, 'CAA94089')
        self.assertEqual(3, len(pr.other_accessions))
        self.assertEqual(set([('embl', 'CAA94089'), ('gi', 'gi|3255965'), ('embl', 'CAA94089.1')]), set(pr.other_accessions))

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
        mock_pm = mock.Mock(spec=Proteome.ProteinManager, uniprotShortcut=False)
        mock_pm.get_protein_sequence.return_value = None

        try:
            entrez_tools.get_protein_information(mock_pm, 'P50914')
        except entrez_tools.EntrezError, e:
            self.assertEqual('P50914', e.acc)
        else:
            self.fail("Expected EntrezError")

    def test_get_protein_information(self):
        pm = Proteome.ProteinManager(settings.adminEmail, uniprotShortcut=False)
        
        pr = entrez_tools.get_protein_information(pm, 'A6NC57')
        
        self.assertEqual('Ankyrin repeat domain-containing protein 62', pr.name)
        self.assertEqual('ANKRD62', pr.gene)
        self.assertEqual('ANR62_HUMAN', pr.locus)
       
        taxons = [u'Eukaryota', u'Metazoa', u'Chordata', u'Craniata',
                u'Vertebrata', u'Euteleostomi', u'Mammalia', u'Eutheria',
                u'Euarchontoglires', u'Primates', u'Haplorrhini',
                u'Catarrhini', u'Hominidae', u'Homo', u'Homo sapiens']

        self.assertEqual('Homo sapiens', pr.species)
        self.assertEqual([t.lower() for t in taxons], pr.taxonomy)
        self.assertEqual('MEVRGSFLAACRRRMATWRKNRDKDGFSNPGYRVRQKDLGMIHKAAIAGDVNKVMESILLRLNDLNDRDKKNRTALLLACAHGRPGVVADLVARKCQLNLTDSENRTALIKAVQCQEEVCASILLEHGANPNVRDMYGNTALHYAIDNENISMARKLLAYGADIEARSQDGHTSLLLAVNRKKEQMVAFLLKKKPDLTAIDNFGRTALILAARNGSTSVVYQLLQHNIDVFCQDISGWTAEDYAVASKFQAIRGMISEYKANKRCKSLQNSNSEQDLEMTSEGEQERLEGCESSQPQVEEKMKKCRNKKMEVSRNVHADDSDNYNDDVDELIHKIKNRKPDNHQSPGKENGEFDRLARKTSNEKSKVKSQIYFTDDLNDISGSSEKTSEDDELPYSDDENFMLLIEQSGMECKDFVSLSKSKNATAACGRSIEDQKCYCERLKVKFQKMKNNISVLQKVLSETDKTKSQSEHQNLQGKKKLCNLRFILQQQEEERIKAEELYEKDIEELKIMEEQYRTQTEVKKQSKLTLKSLEVELKTVRSNSNQNFHTHERERDLWQENHLMRDEIARLRLEIDTIKHQNQETENKYFKDIEIIKENNEDLEKTLKRNEEALTKTITRYSKELNVLMDENTMLNSELQKEKQSMSRLETEMESYRCRLAAALCDHDQRQSSKRDLQLAFQSTVNEWCHLQEDTNSHIQILSQQLSKAESTSSGLETELHYEREALKEKTLHIEHMQGVLSRTQRRLEDIEHMYQNDQPILEKYVRKQQSVEDGLFQLQSQNLLYQQQCNDARKKADNQEKTIINIQVKCEDTVEKLQAECRKLEENNKGLMKECTLLKERQCQYEKEKEEREVVRRQLQREVDDALNKQLLLEAMLEISSERRINLEDEAQSLKKKLGQMRSQVCMKLSMSTVTL', pr.sequence) 
        
        parsed_features = []
        for d in pr.features:
            parsed_features.append( (d.type, d.name, d.start, d.stop, d.source) )
       
        self.assertEqual(4, len(pr.mutations))

        m_test_list = set()
        for m in pr.mutations:
            self.assertEqual('Substitution (single)', m.mutationType)
            self.assertEqual('A6NC57', m.acc_id)
            m_test_list.add(( m.location, m.original, m.mutant ))

        expected_muts = set([(188, 'A', 'S'),
                             (265, 'C', 'R'),
                             (406, 'E', 'K'),
                             (613, 'A', 'T')])

        self.assertEqual(expected_muts, m_test_list)

        exp_features = [ \
                        ('Region', 'Repetitive region', 71, 100, 'ncbi'),
                        ('Region', 'ANK', 99, 224, 'ncbi'),
                        ('Region', 'Repetitive region', 104, 133, 'ncbi'),
                        ('Region', 'Repetitive region', 137, 166, 'ncbi'),
                        ('Region', 'Repetitive region', 170, 199, 'ncbi'),
                        ('Region', 'Repetitive region', 203, 232, 'ncbi'),
                        ('Region', 'Coiled-coil region', 439, 669, 'ncbi'),
                        ('Region', 'Coiled-coil region', 778, 906, 'ncbi'),
                        ]
        self.assertEqual(exp_features, parsed_features)
