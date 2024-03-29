from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker.helpers import upload_helpers
from ptmscout.database import modifications, experiment
from tests.views.mocking import createMockExperiment, createMockProtein,\
    createMockAccession, createMockSpecies,\
    createMockTaxonomy, createMockMeasurement, createMockDomain
from mock import patch
from ptmscout.utils import uploadutils
from ptmscout.config import settings, strings
import os
from ptmscout.database import DBSession

class UploadHelpersTest(IntegrationTestCase):

    def test_find_activation_loop_should_find_loops(self):
        prot = createMockProtein()
        prot.sequence = """
MADNEKLDNQ RLdfgKNKGR DLETMRRQRN EVVVELRKNK RDEHLLKRRN
sldEDICEDS DIDGDYRVQN TSLEAIVQNA SSDNQGIQLS AVQAARKLLS
SDRNPPIdyg IKSGILPILV HCLERDDNPS LQFEAAWALT NIASGTppeT

QAVVQSNAVP LFLRLLHSPH QNVCEQAVWA LGNIIGDGPQ CRDYVISLGV
VKPLLSFISP SIPITFLRNV TWVMVNLCRH KDPPPPMETI QEILPALCVL
IHHTDVNILV DTVdlgSYLT DAGgidIQMV IDSGIVPHLV PLLSHQEVKV

QAVVQSNAVP LFLRLLHSPH QNVCEQAVWA LGNIIGDGPQ CRDYVISLGV
VKdfgSFISP SIPITFLRNV TWVMapeCRH KDPPPPMETI QEILPALCVL
IHHTDVNILV DTVWALSYLT DAGNEQIQMV IDSGIVPHLV PLLSHQEVKV
"""

        prot.sequence = prot.sequence.replace(" ", "").replace("\n","").upper()

        d1 = createMockDomain(prot.id, label='Pkinase')
        d1.start, d1.stop = 10, 60

        d2 = createMockDomain(prot.id, label='Pkinase_Tyr')
        d2.start, d2.stop = 100, 170

        d3 = createMockDomain(prot.id, label='Pkinase')
        d3.start, d3.stop = 260, 330

        d4 = createMockDomain(prot.id, label='Other')
        d4.start, d4.stop = 350, 400

        prot.domains = [ d1, d2, d3, d4 ]

        upload_helpers.find_activation_loops( prot )

        created_regions = [ ( r.type, r.label, r.source, r.start, r.stop ) for r in prot.regions ]

        exp_regions = [( 'Activation Loop', 'Kinase Activation Loop', 'predicted', 13, 53 ),
                       ( 'Activation Loop', 'Possible Kinase Activation Loop', 'predicted', 108, 149 )]

        self.assertEqual( exp_regions, created_regions )


    def test_check_ambiguity_should_find_peptides(self):
        ms = createMockMeasurement(35546, 1)
        ms.peptide = 'DQGsLCTs'
        upload_helpers.check_ambiguity(ms, 'homo sapiens')

        self.assertEqual(5, len(ms.ambiguities))
        set_result = set([ (amb.locus, amb.alt_accession, amb.ms_id) for amb in ms.ambiguities ])

        exp_result = set([('IKKA_HUMAN', 'O15111', ms.id), 
                            ('IKKB_HUMAN', 'O14920', ms.id), 
                            ('IKKB_HUMAN', 'O14920-2', ms.id), 
                            ('IKKB_HUMAN', 'O14920-3', ms.id), 
                            ('IKKB_HUMAN', 'O14920-4', ms.id)])
        self.assertEqual(exp_result, set_result)

    def test_get_taxonomic_lineage_with_taxon_not_found(self):
        exp_lineage = [ 'root node of taxonomy', 'bacteria', 'proteobacteria', 'gammaproteobacteria', 'xanthomonadales', 'xanthomonadaceae',
             'xanthomonas', 'xanthomonas campestris', 'xanthomonas campestris pv. oryzae' ]
        taxon_sp = 'Xanthomonas campestris pv. oryzae'

        taxons = upload_helpers.get_taxonomic_lineage(taxon_sp)

        self.assertEqual(exp_lineage, taxons)

    def test_get_taxonomic_lineage(self):
        human_lineage = [ \
            'root node of taxonomy',
            'eukaryota',
            'metazoa',
            'bilateria',
            'chordata',
            'craniata',
            'vertebrata',
            'gnathostomata',
            'euteleostomi',
            'sarcopterygii',
            'tetrapoda',
            'amniota',
            'mammalia',
            'theria',
            'eutheria',
            'euarchontoglires',
            'primates',
            'haplorrhini',
            'simiiformes',
            'catarrhini',
            'hominoidea',
            'hominidae',
            'homininae',
            'homo',
            'homo sapiens' ]

        lineage = upload_helpers.get_taxonomic_lineage('Homo sapiens')

        self.assertEqual(human_lineage, lineage)

    def test_map_expression_probesets(self):
        prot = createMockProtein()
        prot.species_id = 46
        accessions = [('swissprot', 'O43248'), ('swissprot', 'HXC11_HUMAN'), ('gene_synonym', 'HOX3H'), ('ipi', 'IPI00011610'), ('ipi', 'IPI00011610.1'), ('swissprot', 'O43248.1')]
        
        prot.acc_gene = 'SIRPB1'
        prot.accessions = []
        for acc in accessions:
            prot.accessions.append(createMockAccession(prot.id, value=acc[1], atype=acc[0]))
        
        upload_helpers.map_expression_probesets(prot)

        self.assertEqual(set([6271,6460,16607]), set([probe.id for probe in prot.expression_probes]))
        

    @patch('ptmscout.database.taxonomies.getTaxonByName')
    @patch('ptmscout.database.taxonomies.getSpeciesByName')
    def test_find_or_create_species_should_get_strain_or_isolate(self, patch_getSpecies, patch_getTaxon):
        patch_getSpecies.return_value = None
        taxon = createMockTaxonomy(2759)
        patch_getTaxon.return_value = taxon

        species_name = 'sacchromyces cerevisiae (strain ATC2044)'
        sp = upload_helpers.find_or_create_species(species_name)

        patch_getTaxon.assert_called_once_with('sacchromyces cerevisiae', strain='strain ATC2044')
        self.assertEqual(species_name, sp.name)
        self.assertEqual(taxon.node_id, sp.taxon_id)

    @patch('ptmscout.database.taxonomies.getTaxonByName')
    @patch('ptmscout.database.taxonomies.getSpeciesByName')
    def test_find_or_create_species_should_raise_error_if_no_taxon(self, patch_getSpecies, patch_getTaxon):
        patch_getSpecies.return_value = None
        patch_getTaxon.return_value = None
        
        try:
            upload_helpers.find_or_create_species('some species')
        except uploadutils.ParseError, e:
            self.assertEqual("Species: some species does not match any taxon node", e.message)
        else:
            self.fail("Expected parseerror")

    def test_find_or_create_species_should_query_NCBI_if_no_taxon_found_in_database(self):
        sp = upload_helpers.find_or_create_species('Xanthomonas campestris pv. oryzae')

        self.assertEqual(sp.taxon_id, 314227)
        self.assertEqual(sp.name, 'Xanthomonas campestris pv. oryzae')

    @patch('ptmscout.database.taxonomies.getTaxonByName')
    @patch('ptmscout.database.taxonomies.getSpeciesByName')
    def test_find_or_create_species_should_create_species_if_taxon_node_available(self, patch_getSpecies, patch_getTaxon):
        patch_getSpecies.return_value = None
        taxon = createMockTaxonomy(2759)
        patch_getTaxon.return_value = taxon
        
        sp = upload_helpers.find_or_create_species('some species')
        
        patch_getTaxon.assert_called_once_with('some species', strain=None)
        self.assertEqual('some species', sp.name)
        self.assertEqual(taxon.node_id, sp.taxon_id)

    @patch('ptmscout.database.taxonomies.getSpeciesByName')
    def test_find_or_create_species_should_return_species_if_available(self, patch_getSpecies):
        exp_species = createMockSpecies()
        patch_getSpecies.return_value = exp_species
        
        sp = upload_helpers.find_or_create_species('some species')
        self.assertEqual(exp_species, sp)

    @patch('ptmscout.database.modifications.Peptide')    
    @patch('ptmscout.database.modifications.getPeptideBySite')
    def test_get_peptide_should_create_new_peptide_if_not_found(self, patch_getPep, patch_pep):
        pid = 100
        pep_site = 200
        pep_seq = "ASDFGBSdKEDFASD"
        pep_type = 'D'
                
        patch_getPep.return_value = None
        patch_getPep.side_effect = modifications.NoSuchPeptide(pep_site, pep_type, pid)
        
        pepinst, created =  upload_helpers.get_peptide(pid, pep_site, pep_seq)
        
        self.assertTrue(created)        
        self.assertEqual(pep_seq, pepinst.pep_aligned)
        self.assertEqual(pep_site, pepinst.site_pos)
        self.assertEqual(pep_type, pepinst.site_type)
        self.assertEqual(pid, pepinst.protein_id)
        
    def test_insert_run_data_when_exists_should_modify_existing(self):
        MS_peptide = modifications.MeasuredPeptide()
        MS_peptide.query_accession = 'Q01234'
        MS_peptide.experiment_id = 1
        MS_peptide.protein_id = 35546
        MS_peptide.peptide = 'ABCDEFG'
        
        data1 = experiment.ExperimentData()
        data1.run = 'run1'
        data1.label = '0'
        data1.type = 'data'
        data1.priority = 4
        data1.value = 9

        data2 = experiment.ExperimentData()
        data2.run = 'run1'
        data2.label = '20'
        data2.type = 'data'
        data2.priority = 3
        data2.value = 10



        MS_peptide.data.extend([data1, data2])

        DBSession.add(MS_peptide)
        DBSession.flush()
        
        series_header = [('data', '0'),('data', '5'),('data', '20'), ('stddev', '5'), ('stddev', '20')]
        series = [0,4,1,3,2]
        
        upload_helpers.insert_run_data(MS_peptide, 1, 'time(min)', series_header, "run1", series)
        
        DBSession.flush()
        
        result = DBSession.query(experiment.ExperimentData).filter_by(MS_id=MS_peptide.id).all()
        result = sorted(result, key=lambda item: item.priority)
        
        self.assertEqual(5, len(result))
        
        self.assertEqual("run1", result[0].run)
        self.assertEqual("data", result[0].type)
        self.assertEqual("time(min)", result[0].units)
        self.assertEqual('0', result[0].label)
        self.assertEqual(0, result[0].value)
        
        self.assertEqual("run1", result[1].run)
        self.assertEqual("data", result[1].type)
        self.assertEqual("time(min)", result[1].units)
        self.assertEqual('5', result[1].label)
        self.assertEqual(4, result[1].value)

        self.assertEqual("run1", result[2].run)
        self.assertEqual("data", result[2].type)
        self.assertEqual("time(min)", result[2].units)
        self.assertEqual('20', result[2].label)
        self.assertEqual(1, result[2].value)

        self.assertEqual("run1", result[3].run)
        self.assertEqual("stddev", result[3].type)
        self.assertEqual("time(min)", result[3].units)
        self.assertEqual('5', result[3].label)
        self.assertEqual(3, result[3].value)

        self.assertEqual("run1", result[4].run)
        self.assertEqual("stddev", result[4].type)
        self.assertEqual("time(min)", result[4].units)
        self.assertEqual('20', result[4].label)
        self.assertEqual(2, result[4].value)
        
 

    def test_insert_run_data_should_create_data_records_when_timeseries(self):
        MS_peptide = modifications.MeasuredPeptide()
        MS_peptide.experiment_id = 1
        MS_peptide.protein_id = 35546
        MS_peptide.query_accession = 'Q01234'
        MS_peptide.peptide = 'ABCDEFG'
        

        DBSession.add(MS_peptide)
        DBSession.flush()
        
        series_header = [('data', '0'),('data', '5'),('data', '20'), ('stddev', '5'), ('stddev', '20')]
        series = [0,None,1,3,2]
        
        upload_helpers.insert_run_data(MS_peptide, 1, 'time(min)', series_header, "run1", series)
        
        DBSession.flush()
        
        result = DBSession.query(experiment.ExperimentData).filter_by(MS_id=MS_peptide.id).all()
        result = sorted(result, key=lambda item: item.priority)
        
        self.assertEqual(5, len(result))
        
        self.assertEqual("run1", result[0].run)
        self.assertEqual("data", result[0].type)
        self.assertEqual("time(min)", result[0].units)
        self.assertEqual('0', result[0].label)
        self.assertEqual(0, result[0].value)
        
        self.assertEqual("run1", result[1].run)
        self.assertEqual("data", result[1].type)
        self.assertEqual("time(min)", result[1].units)
        self.assertEqual('5', result[1].label)
        self.assertTrue(result[1].value == None)

        self.assertEqual("run1", result[2].run)
        self.assertEqual("data", result[2].type)
        self.assertEqual("time(min)", result[2].units)
        self.assertEqual('20', result[2].label)
        self.assertEqual(1, result[2].value)

        self.assertEqual("run1", result[3].run)
        self.assertEqual("stddev", result[3].type)
        self.assertEqual("time(min)", result[3].units)
        self.assertEqual('5', result[3].label)
        self.assertEqual(3, result[3].value)

        self.assertEqual("run1", result[4].run)
        self.assertEqual("stddev", result[4].type)
        self.assertEqual("time(min)", result[4].units)
        self.assertEqual('20', result[4].label)
        self.assertEqual(2, result[4].value)
        
   
    def test_parse_modifications_should_find_peptides_at_start(self):
        prot_seq = \
"""MKPTRSQLDSDFSQQDTPCLIIEDSQPESQVLEDDSGSHFSMLSRHLPNLQTHKENPVLD
VVSNPEQTAGEERGDGNSGFNEHLKENKVADPVDSSNLDTCGSISQVIEQLPQPNRTSSV
LGMSVESAPAVEEEKGEELEQKEKEKEEDTSGNTTHSLGAEDTASSQLGFGVLELSQSQD
VEENTVPYEVDKEQLQSVTTNSGYTRLSDVDANTAIKHEEQSNEDIPIAEQSSKDIPVTA
QPSKDVHVVKEQNPPPARSEDMPFSPKASVAAMEAKEQLSAQELMESGLQIQKSPEPEVL
STQEDLFDQSNKTVSSDGCSTPSREEGGCSLASTPATTLHLLQLSGQRSLVQDSLSTNSS
DLVAPSPDAFRSTPFIVPSSPTEQEGRQDKPMDTSVLSEEGGEPFQKKLQSGEPVELENP
PLLPESTVSPQASTPISQSTPVFPPGSLPIPSQPQFSHDIFIPSPSLEEQSNDGKKDGDM
HSSSLTVECSKTSEIEPKNSPEDLGLSLTGDSCKLMLSTSEYSQSPKMESLSSHRIDEDG
ENTQIEDTEPMSPVLNSKFVPAENDSILMNPAQDGEVQLSQNDDKTKGDDTDTRDDISIL
"""
        prot_seq = prot_seq.replace("\n", "")
        pep_seq = "  MkPTrSQLDSDF"
        
        taxonomy = set(['eukaryota', 'chordata', 'mammalia', 'homo'])
        mod_types, aligned_peps = upload_helpers.parse_modifications(prot_seq, pep_seq, "METHYLATION", taxonomy)
        
        self.assertEqual('N6-methyllysine', mod_types[0].name)
        self.assertEqual('Omega-N-methylarginine', mod_types[1].name)
       
        self.assertEqual((2, "      MkPTRSQLD", 'K'), aligned_peps[0])
        self.assertEqual((5, "   MKPTrSQLDSDF", 'R'), aligned_peps[1])

    def test_protein_sequence_match_ambiguous_protein(self):
        prot_seq = \
"""
TSDKLASRSKLPDGPTGSSEEEEEFLEIPPFNKQYTESQLRAGAGYILEDFNEAQCNTAY
QCLLIADQHCRTRKYFLCLASGIPCVSHVWVHDSCHANQLQNYRNYLLPAGYSLEEQRIL
QCLLIADQHCRTRKYFLCLASGIPCVSHVWVHDSCHANQLQNYRNYLLPAGYSLEEQRIL
DWQPRENPFQNLKVLLVSDQQQNFLELWSEILMTGGAASVKQHHSSAHNKDIALGVFDVV
VTDPSCPASVLKCAEALQLPVVSQEWVIQCLIVGERIGFKQHPKYKHD"""        
        prot_seq = prot_seq.replace("\n", "")
        pep_seq = "CHANQLQNYRN"
        
        taxonomy = set(['eukaryota', 'chordata', 'mammalia', 'homo'])
        
        try:
            upload_helpers.parse_modifications(prot_seq, pep_seq, "METHYLATION", taxonomy)
        except uploadutils.ParseError, e:
            self.assertEqual(strings.experiment_upload_warning_peptide_ambiguous_location_in_protein_sequence, e.msg)
        else:
            self.fail("Expected parse error")
        

    def test_protein_sequence_match_end_of_protein(self):
        prot_seq = \
"""
TSDKLASRSKLPDGPTGSSEEEEEFLEIPPFNKQYTESQLRAGAGYILEDFNEAQCNTAY
QCLLIADQHCRTRKYFLCLASGIPCVSHVWVHDSCHANQLQNYRNYLLPAGYSLEEQRIL
DWQPRENPFQNLKVLLVSDQQQNFLELWSEILMTGGAASVKQHHSSAHNKDIALGVFDVV
VTDPSCPASVLKCAEALQLPVVSQEWVIQCLIVGERIGFKQHPKYKHD"""        
        prot_seq = prot_seq.replace("\n", "")
        pep_seq = "ERIGFKQHPkYkHD      "
        
        taxonomy = set(['eukaryota', 'chordata', 'mammalia', 'homo'])
        mod_types, aligned_peps = upload_helpers.parse_modifications(prot_seq, pep_seq, "METHYLATION", taxonomy)
        
        self.assertEqual('N6-methyllysine', mod_types[0].name)
        self.assertEqual('N6-methyllysine', mod_types[1].name)
       
        self.assertEqual((224, "IGFKQHPkYKHD   ", 'K'), aligned_peps[0])
        self.assertEqual((226, "FKQHPKYkHD     ", 'K'), aligned_peps[1])


    

    def test_get_aligned_peptide_sequences_should_produce_correct_seqs(self):
        pep_seq = "mSaDFJTkLJAWERpOID"
        pep_up = pep_seq.upper()
        prot_seq = "MSADFJTKLJAWERPOIDFK"
        
        mod_sites = [ i for i in xrange(0, len(pep_seq)) if pep_seq[i] != pep_up[i] ]
        index = 0
        
        aligned_peptides = upload_helpers.get_aligned_peptide_sequences(mod_sites, index, pep_seq, prot_seq)
        
        self.assertEqual((1,  "       mSADFJTK", 'M'), aligned_peptides[0])
        self.assertEqual((3,  "     MSaDFJTKLJ", 'A'), aligned_peptides[1])
        self.assertEqual((8,  "MSADFJTkLJAWERP", 'K'), aligned_peptides[2])
        self.assertEqual((15, "KLJAWERpOIDFK  ", 'P'), aligned_peptides[3])



    def test_parse_modifications_should_produce_correct_seqs(self):
        prot_seq = \
"""MDPTGSQLDSDFSQQDTPCLIIEDSQPESQVLEDDSGSHFSMLSRHLPNLQTHKENPVLD
VVSNPEQTAGEERGDGNSGFNEHLKENKVADPVDSSNLDTCGSISQVIEQLPQPNRTSSV
LGMSVESAPAVEEEKGEELEQKEKEKEEDTSGNTTHSLGAEDTASSQLGFGVLELSQSQD
VEENTVPYEVDKEQLQSVTTNSGYTRLSDVDANTAIKHEEQSNEDIPIAEQSSKDIPVTA
QPSKDVHVVKEQNPPPARSEDMPFSPKASVAAMEAKEQLSAQELMESGLQIQKSPEPEVL
STQEDLFDQSNKTVSSDGCSTPSREEGGCSLASTPATTLHLLQLSGQRSLVQDSLSTNSS
DLVAPSPDAFRSTPFIVPSSPTEQEGRQDKPMDTSVLSEEGGEPFQKKLQSGEPVELENP
PLLPESTVSPQASTPISQSTPVFPPGSLPIPSQPQFSHDIFIPSPSLEEQSNDGKKDGDM
HSSSLTVECSKTSEIEPKNSPEDLGLSLTGDSCKLMLSTSEYSQSPKMESLSSHRIDEDG
ENTQIEDTEPMSPVLNSKFVPAENDSILMNPAQDGEVQLSQNDDKTKGDDTDTRDDISIL
ATGCKGREETVAEDVCIDLTCDSGSQAVPSPATRSEALSSVLDQEEAMEIKEHHPEEGSS
GSEVEEIPETPCESQGEELKEENMESVPLHLSLTETQSQGLCLQKEMPKKECSEAMEVET
SVISIDSPQKLAILDQELEHKEQEAWEEATSEDSSVVIVDVKEPSPRVDVSCEPLEGVEK
CSDSQSWEDIAPEIEPCAENRLDTKEEKSVEYEGDLKSGTAETEPVEQDSSQPSLPLVRA
DDPLRLDQELQQPQTQEKTSNSLTEDSKMANAKQLSSDAEAQKLGKPSAHASQSFCESSS
ETPFHFTLPKEGDIIPPLTGATPPLIGHLKLEPKRHSTPIGISNYPESTIATSDVMSESM
VETHDPILGSGKGDSGAAPDVDDKLCLRMKLVSPETEASEESLQFNLEKPATGERKNGST
AVAESVASPQKTMSVLSCICEARQENEARSEDPPTTPIRGNLLHFPSSQGEEEKEKLEGD
HTIRQSQQPMKPISPVKDPVSPASQKMVIQGPSSPQGEAMVTDVLEDQKEGRSTNKENPS
KALIERPSQNNIGIQTMECSLRVPETVSAATQTIKNVCEQGTSTVDQNFGKQDATVQTER
GSGEKPVSAPGDDTESLHSQGEEEFDMPQPPHGHVLHRHMRTIREVRTLVTRVITDVYYV
DGTEVERKVTEETEEPIVECQECETEVSPSQTGGSSGDLGDISSFSSKASSLHRTSSGTS
LSAMHSSGSSGKGAGPLRGKTSGTEPADFALPSSRGGPGKLSPRKGVSQTGTPVCEEDGD
AGLGIRQGGKAPVTPRGRGRRGRPPSRTTGTRETAVPGPLGIEDISPNLSPDDKSFSRVV
PRVPDSTRRTDVGAGALRRSDSPEIPFQAAAGPSDGLDASSPGNSFVGLRVVAKWSSNGY
FYSGKITRDVGAGKYKLLFDDGYECDVLGKDILLCDPIPLDTEVTALSEDEYFSAGVVKG
HRKESGELYYSIEKEGQRKWYKRMAVILSLEQGNRLREQYGLGPYEAVTPLTKAADISLD
NLVEGKRKRRSNVSSPATPTASSSSSTTPTRKITESPRASMGVLSGKRKLITSEEERSPA
KRGRKSATVKPGAVGAGEFVSPCESGDNTGEPSALEEQRGPLPLNKTLFLGYAFLLTMAT
TSDKLASRSKLPDGPTGSSEEEEEFLEIPPFNKQYTESQLRAGAGYILEDFNEAQCNTAY
QCLLIADQHCRTRKYFLCLASGIPCVSHVWVHDSCHANQLQNYRNYLLPAGYSLEEQRIL
DWQPRENPFQNLKVLLVSDQQQNFLELWSEILMTGGAASVKQHHSSAHNKDIALGVFDVV
VTDPSCPASVLKCAEALQLPVVSQEWVIQCLIVGERIGFKQHPKYKHDYVSH"""        
        prot_seq = prot_seq.replace("\n", "")
        pep_seq = "GkAPVTPrGrGrrGr"
        
        taxonomy = set(['eukaryota', 'chordata', 'mammalia', 'homo'])
        mod_types, aligned_peps = upload_helpers.parse_modifications(prot_seq, pep_seq, "METHYLATION", taxonomy)
        
        self.assertEqual('N6-methyllysine', mod_types[0].name)
        self.assertEqual('Omega-N-methylarginine', mod_types[1].name)
        self.assertEqual('Omega-N-methylarginine', mod_types[2].name)
        self.assertEqual('Omega-N-methylarginine', mod_types[3].name)
        self.assertEqual('Omega-N-methylarginine', mod_types[4].name)
        self.assertEqual('Omega-N-methylarginine', mod_types[5].name)
        
        self.assertEqual((1390, "LGIRQGGkAPVTPRG", 'K'), aligned_peps[0])
        self.assertEqual((1396, "GKAPVTPrGRGRRGR", 'R'), aligned_peps[1])
        self.assertEqual((1398, "APVTPRGrGRRGRPP", 'R'), aligned_peps[2])
        self.assertEqual((1400, "VTPRGRGrRGRPPSR", 'R'), aligned_peps[3])
        self.assertEqual((1401, "TPRGRGRrGRPPSRT", 'R'), aligned_peps[4])
        self.assertEqual((1403, "RGRGRRGrPPSRTTG", 'R'), aligned_peps[5])
    
    def test_get_aligned_peptide_sequences_should_produce_correct_seqs(self):
        pep_seq = "mSaDFJTkLJAWERpOID"
        pep_up = pep_seq.upper()
        prot_seq = "MSADFJTKLJAWERPOIDFK"
        
        mod_sites = [ i for i in xrange(0, len(pep_seq)) if pep_seq[i] != pep_up[i] ]
        index = 0
        
        aligned_peptides = upload_helpers.get_aligned_peptide_sequences(mod_sites, index, pep_seq, prot_seq)
        
        self.assertEqual((1,  "       mSADFJTK", 'M'), aligned_peptides[0])
        self.assertEqual((3,  "     MSaDFJTKLJ", 'A'), aligned_peptides[1])
        self.assertEqual((8,  "MSADFJTkLJAWERP", 'K'), aligned_peptides[2])
        self.assertEqual((15, "KLJAWERpOIDFK  ", 'P'), aligned_peptides[3])


    def test_create_chunked_tasks_preserve_groups_should_build_correct_tasks(self):
        task_args = ['A53D56', 'F435D6', 'ALB432', 'Q134A5', 'Q134A5-3', 'Q134A5-5', 'P54A03', 'E45G76', 'E45G76-2']
        
        tasks = upload_helpers.create_chunked_tasks_preserve_groups(sorted(task_args), 4)
        
        self.assertEqual(['A53D56', 'ALB432', 'E45G76', 'E45G76-2'], tasks[0])
        self.assertEqual(['F435D6', 'P54A03'], tasks[1])
        self.assertEqual(['Q134A5', 'Q134A5-3', 'Q134A5-5'], tasks[2])

    def test_create_chunked_tasks_should_build_correct_tasks(self):
        task_args = ['A53D56', 'F435D6', 'ALB432', 'Q134A5', 'Q134A5-3', 'Q134A5-5', 'P54A03', 'E45G76', 'E45G76-2']
        
        tasks = upload_helpers.create_chunked_tasks(sorted(task_args), 4)
        
        self.assertEqual(['A53D56', 'ALB432', 'E45G76', 'E45G76-2'], tasks[0])
        self.assertEqual(['F435D6', 'P54A03', 'Q134A5', 'Q134A5-3'], tasks[1])
        self.assertEqual(['Q134A5-5'], tasks[2])

    def test_store_last_stage_result_should_save_result_to_file(self):
        exp = createMockExperiment()
        exp.id = 1002304052435
        exp.loading_stage = 'GO terms'
        result = {"some":"pyobject", "which":"canbe", "pickled":"."}

        upload_helpers.store_stage_input(exp.id, exp.loading_stage, result)

        result_path = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'e%d') % (exp.id)
        result_file = 'GO_terms.input'

        loaded_result = upload_helpers.get_stage_input(exp.id, exp.loading_stage)

        self.assertEqual(result, loaded_result)

        os.remove(os.path.join(result_path, result_file))
        os.removedirs(result_path)
