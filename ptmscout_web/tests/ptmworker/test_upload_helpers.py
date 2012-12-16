from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker import upload_helpers
from ptmscout.database import modifications, experiment
from tests.views.mocking import createMockExperiment, createMockProtein,\
    createMockProbe, createMockAccession, createMockSpecies, createMockTaxonomy
from mock import patch
from ptmscout.utils import uploadutils

class PTMWorkerUploadHelpersTestCase(IntegrationTestCase):

    @patch('ptmscout.database.gene_expression.getExpressionProbeSetsForProtein')
    def test_map_expression_probesets(self, patch_getProbes):
        probe = createMockProbe()
        prot = createMockProtein()
        
        prot.acc_gene = 'TNK2'
        prot.accessions.append(createMockAccession(prot.id, value='ACK1_HUMAN', type='uniprot'))
        
        patch_getProbes.return_value = [probe]
        
        upload_helpers.map_expression_probesets(prot)
        
        patch_getProbes.assert_called_once_with(['ACK1_HUMAN', 'TNK2'], prot.species_id)
        
        self.assertEqual([probe], prot.expression_probes)

    @patch('ptmscout.database.taxonomies.getTaxonByName')
    @patch('ptmscout.database.taxonomies.getSpeciesByName')
    def test_find_or_create_species_should_raise_error_if_no_taxon(self, patch_getSpecies, patch_getTaxon):
        patch_getSpecies.return_value = None
        taxon = createMockTaxonomy()
        patch_getTaxon.return_value = None
        
        try:
            upload_helpers.find_or_create_species('some species')
        except uploadutils.ParseError, e:
            self.assertEqual("Species: some species does not match any taxon node", e.message)
        else:
            self.fail("Expected parseerror")
 
    @patch('ptmscout.database.taxonomies.getTaxonByName')
    @patch('ptmscout.database.taxonomies.getSpeciesByName')
    def test_find_or_create_species_should_create_species_if_taxon_node_available(self, patch_getSpecies, patch_getTaxon):
        patch_getSpecies.return_value = None
        taxon = createMockTaxonomy()
        patch_getTaxon.return_value = taxon
        
        sp = upload_helpers.find_or_create_species('some species')
        
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
        pepinst.save.assert_called_once_with()
        
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_mark_experiment_should_get_experiment_and_mark(self, patch_getExp):
        exp = createMockExperiment()
        exp.status = 'loading'
        patch_getExp.return_value = exp
        
        v = upload_helpers.mark_experiment(exp.id, 'loaded')
        
        patch_getExp.assert_called_once_with(exp.id, check_ready=False, secure=False)
        
        self.assertEqual('loaded', exp.status)
        exp.saveExperiment.assert_called_once_with()
        self.assertEqual(exp, v)

    def test_insert_run_data_when_exists_should_modify_existing(self):
        from ptmscout.database import DBSession
        MS_peptide = modifications.MeasuredPeptide()
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
        from ptmscout.database import DBSession
        MS_peptide = modifications.MeasuredPeptide()
        MS_peptide.experiment_id = 1
        MS_peptide.protein_id = 35546
        MS_peptide.peptide = 'ABCDEFG'
        

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
        class DummyTask(object):
            def __init__(self, args=None):
                self.args = args

            def s(self, *args):
                return DummyTask(args=args)

        task_args = ['A53D56', 'F435D6', 'ALB432', 'Q134A5', 'Q134A5-3', 'Q134A5-5', 'P54A03', 'E45G76', 'E45G76-2']
        
        tasks = upload_helpers.create_chunked_tasks_preserve_groups(DummyTask(), sorted(task_args), 4)
        
        self.assertEqual((['A53D56', 'ALB432', 'E45G76', 'E45G76-2'],), tasks[0].args)
        self.assertEqual((['F435D6', 'P54A03'],), tasks[1].args)
        self.assertEqual((['Q134A5', 'Q134A5-3', 'Q134A5-5'],), tasks[2].args)

    def test_create_chunked_tasks_should_build_correct_tasks(self):
        class DummyTask(object):
            def __init__(self, args=None):
                self.args = args

            def s(self, *args):
                return DummyTask(args=args)

        task_args = ['A53D56', 'F435D6', 'ALB432', 'Q134A5', 'Q134A5-3', 'Q134A5-5', 'P54A03', 'E45G76', 'E45G76-2']
        
        tasks = upload_helpers.create_chunked_tasks(DummyTask(), sorted(task_args), 4)
        
        self.assertEqual((['A53D56', 'ALB432', 'E45G76', 'E45G76-2'],), tasks[0].args)
        self.assertEqual((['F435D6', 'P54A03', 'Q134A5', 'Q134A5-3'],), tasks[1].args)
        self.assertEqual((['Q134A5-5'],), tasks[2].args)
