from tests.PTMScoutTestCase import IntegrationTestCase
from ptmscout.config import settings
from Bio import Entrez
import os
from ptmworker import entrez_tools
import json
import mock
from geeneus import Proteome

class EntrezQueryTestCase(IntegrationTestCase):

    def test_get_pubmed_record(self):
        record = entrez_tools.get_pubmed_record_by_id(12230038)
        
        self.assertEqual(settings.adminEmail, Entrez.email)
        
        self.assertEqual('The Bio* toolkits--a brief overview.', record['TI'])
        self.assertEqual('296-302', record['PG'])
        self.assertEqual('2002 Sep', record['DP'])


    def test_parse_gene_pfam_domains_should_return_parsed_features(self):
        path = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'json_proteinfeaturexml')
        proteinfeature_xml = json.loads(open(path, 'rb').read())
        
        protein_domains = entrez_tools.parse_proteinxml_pfam_domains(proteinfeature_xml)
        
        self.assertEqual(1, len(protein_domains))
        
        self.assertEqual(1, protein_domains[0].significant)
        self.assertEqual(0, protein_domains[0].release)
        self.assertEqual(-1, protein_domains[0].p_value)
        self.assertEqual('Ribosomal_L14e', protein_domains[0].label)
        self.assertEqual(47, protein_domains[0].start)
        self.assertEqual(119, protein_domains[0].stop)
        self.assertEqual('Domain', protein_domains[0].class_)
    
    def test_parse_gene_name_should_return_name(self):
        path = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'json_proteinfeaturexml')
        proteinfeature_xml = json.loads(open(path, 'rb').read())
        gene_name = entrez_tools.get_gene_name(proteinfeature_xml)
        
        self.assertEqual("RPL14", gene_name)
        
    def test_get_protein_information_should_raise_entrez_error(self):
        mock_pm = mock.Mock(spec=Proteome.ProteinManager)
        mock_pm.get_protein_sequence.return_value = ""

        try:
            entrez_tools.get_protein_information(mock_pm, 'P50914')
        except entrez_tools.EntrezError, e:
            self.assertEqual('P50914', e.acc)
        else:
            self.fail("Expected EntrezError")

    def test_get_protein_information_should_return_correct_info(self):
        path = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'json_proteinxml')
        protein_xml = json.loads(open(path, 'rb').read())
        
        pm = mock.Mock(spec=Proteome.ProteinManager)
        exp_seq = 'mvfrrfvevgrvayvsfgphagklvaivdvidqnralvdgpctqvrrqampfkcmqltdfilkfphsahqkyvrqawqkadintkwaatrwakkiearerkakmtdfdrfkvmkakkmrnriiknevkklqkaallkaspkkapgtkgtaaaaaaaaaakvpakkitaaskkapaqkvpaqkatgqkaapapkaqkgqkapaqkapapkasgkka'
        pm.get_protein_sequence.return_value = exp_seq
        pm.get_raw_xml.return_value = protein_xml
        
        prot_name, gene, taxonomy, species, prot_accessions, prot_domains, seq = entrez_tools.get_protein_information(pm, 'P50914')
        
        self.assertEqual('RecName: Full=60S ribosomal protein L14; AltName: Full=CAG-ISL 7', prot_name)
        self.assertEqual('RPL14', gene)
        
        expected_taxonomy = 'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo'.split("; ")
        self.assertEqual(set([ tax.lower() for tax in expected_taxonomy ]), taxonomy)
        
        self.assertEqual('Homo sapiens', species)
        self.assertEqual(1, len(prot_domains))
        
        self.assertEqual(exp_seq.upper(), seq)
        
        expected_accessions = set([('swissprot', 'RL14_HUMAN'), ('swissprot', 'P50914'), ('gi','gi|212276521')])
        self.assertEqual(expected_accessions, prot_accessions)
        
