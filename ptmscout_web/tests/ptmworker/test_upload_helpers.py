from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker import upload_helpers
from ptmscout.database import modifications, experiment
from tests.views.mocking import createMockExperiment
from mock import patch

class PTMWorkerUploadHelpersTestCase(IntegrationTestCase):

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
            