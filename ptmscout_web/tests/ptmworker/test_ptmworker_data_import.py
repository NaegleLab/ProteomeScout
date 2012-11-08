# Skip
from tests.PTMScoutTestCase import IntegrationTestCase
from mock import patch, call
from ptmworker import tasks
import os
from ptmscout.database import experiment, modifications, protein
from ptmscout.config import settings
from tests.views.mocking import createMockExperiment

class PTMWorkDataImportTestCase(IntegrationTestCase):
    
    def test_finalize_import_should_set_experiment_to_loaded_status(self):
        exp = createMockExperiment()
        exp.status = 'loading'
        
        res = tasks.finalize_import.apply_async((exp,))
        res.get()
        assert res.successful()
        
        exp.status = 'loaded'
        exp.saveExperiment.assert_called_once_with()
    
    def test_load_peptide_should_fail_to_align_if_sequence_is_missing(self):
        prot = protein.getProteinById(35546)
        
        # MQPEE GTGWL LELLS EVQLQ QYFLR LRDDL NVTRL--SHFEY VKNED LEKIG M--GRPG QRRLW EAVKRRKALCKRKSWMSKVFSGKRLEAEFPPHHSQSTFRKTSPAPGGPAGEGPLQSLTCLIGEKDLRLLEKLGDGSFGVVRRGEWDAPSGKTVSVAVKCLKPDVLSQPEAMDDFIREVNAMHSLDHRNLIRLYGVVLTPPMKMVTELAPLGSLLDRLRKHQGHFLLGTLSRYAVQVAEGMGYLESKRFIHRDLAARNLLLATRDLVKIGDFGLMRALPQNDDHYVMQEHRKVPFAWCAPESLKTRTFSHASDTWMFGVTLWEMFTYGQEPWIGLNGSQILHKIDKEGERLPRPEDCPQDIYNVMVQCWAHKPEDRPTFVALRDFLLEAQPTDMRALQDFEEPDKLHIQMNDVITVIEGRAENYWWRGQNTRTLCVGPFPRNVVTSVAGLSAQDISQPLQNSFIHTGHGDSDPRHCWGFPDRIDELYLGNPMDPPDLLSVELSTSRPPQHLGGVKKPTYDPVSEDQDPLSSDFKRLGLRKPGLPRGLWLAKPSARVPGTKASRGSGAEVTLIDFGEEPVVPALRPCPPSLAQLAMDACSLLDETPPQSPTRALPRPLHPTPVVDWDARPLPPPPAYDDVAQDEDDFEICSINSTLVGAGVPAGPSQGQTNYAFVPEQARPPPPLEDNLFLPPQGGGKPPSSAQTAEIFQALQQECMRQLQAPGSPAPSPSPGGDDKPQVPPRVPIPPRPTRPHVQLSPAPPGEEETSQWPGPASPPRVPPREPLSPQGSRTPSPLVPPGSSPLPPRLSSSPGKTMPTTQSFASDPKYATPQVIQAPGAGGPCILPIVRDGKKVSSTHYYLLPERPSYLERYQRFLREAQSPEEPTPLPVPLLLPPPSTPAPAAPTATVRPMPQAALDPKANFSTNNSNPGARPPPPRATARLPQRGCPGDGPEAGRPADKIQMAMVHGVTTEECQAALQCHGWSVQRAAQYLKVEQLFGLGLRPRGECHKVLEMFDWNLEQAGCHLLGSWGPAHHKR
        m = {'gi|8922075': prot}
        t = tasks.load_peptide.apply_async((m, 1, 'gi|8922075', 'AAAAAAAAAAAaAaAAAQRRL'))
        try:
            t.get()
        except:
            assert t.failed()
        
    
    def test_load_peptide_should_align_peptide_sequence_and_create_record(self):
        from ptmscout.database import DBSession
        
        prot = protein.getProteinById(35546)
        
        ppep = modifications.Phosphopep()
        ppep.pep_tryps = 'LEKiGMGRPGQRRLW'
        ppep.pep_aligned = 'KNEDLEKiGMGRPGQ'
        ppep.protein_id = prot.id
        ppep.pfam_site = "~~~"
        ppep.site_pos = 49
        ppep.site_type = 'I'
        
        DBSession.add(ppep)
        DBSession.flush()
        
        # MQPEE GTGWL LELLS EVQLQ QYFLR LRDDL NVTRL--SHFEY VKNED LEKIG M--GRPG QRRLW EAVKRRKALCKRKSWMSKVFSGKRLEAEFPPHHSQSTFRKTSPAPGGPAGEGPLQSLTCLIGEKDLRLLEKLGDGSFGVVRRGEWDAPSGKTVSVAVKCLKPDVLSQPEAMDDFIREVNAMHSLDHRNLIRLYGVVLTPPMKMVTELAPLGSLLDRLRKHQGHFLLGTLSRYAVQVAEGMGYLESKRFIHRDLAARNLLLATRDLVKIGDFGLMRALPQNDDHYVMQEHRKVPFAWCAPESLKTRTFSHASDTWMFGVTLWEMFTYGQEPWIGLNGSQILHKIDKEGERLPRPEDCPQDIYNVMVQCWAHKPEDRPTFVALRDFLLEAQPTDMRALQDFEEPDKLHIQMNDVITVIEGRAENYWWRGQNTRTLCVGPFPRNVVTSVAGLSAQDISQPLQNSFIHTGHGDSDPRHCWGFPDRIDELYLGNPMDPPDLLSVELSTSRPPQHLGGVKKPTYDPVSEDQDPLSSDFKRLGLRKPGLPRGLWLAKPSARVPGTKASRGSGAEVTLIDFGEEPVVPALRPCPPSLAQLAMDACSLLDETPPQSPTRALPRPLHPTPVVDWDARPLPPPPAYDDVAQDEDDFEICSINSTLVGAGVPAGPSQGQTNYAFVPEQARPPPPLEDNLFLPPQGGGKPPSSAQTAEIFQALQQECMRQLQAPGSPAPSPSPGGDDKPQVPPRVPIPPRPTRPHVQLSPAPPGEEETSQWPGPASPPRVPPREPLSPQGSRTPSPLVPPGSSPLPPRLSSSPGKTMPTTQSFASDPKYATPQVIQAPGAGGPCILPIVRDGKKVSSTHYYLLPERPSYLERYQRFLREAQSPEEPTPLPVPLLLPPPSTPAPAAPTATVRPMPQAALDPKANFSTNNSNPGARPPPPRATARLPQRGCPGDGPEAGRPADKIQMAMVHGVTTEECQAALQCHGWSVQRAAQYLKVEQLFGLGLRPRGECHKVLEMFDWNLEQAGCHLLGSWGPAHHKR
        m = {'gi|8922075': prot}
        t = tasks.load_peptide.apply_async((m, 1, 'gi|8922075', 'SHFEYVKNEdLEKiGM'))
        try:
            MSpeptide = t.get()
        except:
            print t.traceback
        
        assert t.successful()
        
        self.assertEqual('SHFEYVKNEdLEKiGM', MSpeptide.phosphopep)
        self.assertEqual(35546, MSpeptide.protein_id)
        self.assertEqual(1, MSpeptide.experiment_id)
        
        DBSession.flush()
        MSpep = DBSession.query(modifications.MeasuredPeptide).filter_by(id=MSpeptide.id).first()
        
        pep1, pep2 = MSpep.phosphopeps
        
        self.assertEqual('SHFEYVKNEdLEKIGM', pep1.pep_tryps)
        self.assertEqual('FEYVKNEdLEKIGMG', pep1.pep_aligned)
        self.assertEqual(45, pep1.site_pos)
        self.assertEqual('D', pep1.site_type)

        self.assertEqual('LEKiGMGRPGQRRLW', pep2.pep_tryps)
        self.assertEqual('KNEDLEKiGMGRPGQ', pep2.pep_aligned)
        self.assertEqual(49, pep2.site_pos)
        self.assertEqual('I', pep2.site_type)
    
    def test_insert_run_data_should_create_data_records_when_timeseries(self):
        from ptmscout.database import DBSession
        MS_peptide = modifications.MeasuredPeptide()
        MS_peptide.experiment_id = 1
        MS_peptide.protein_id = 35546
        MS_peptide.phosphopep = 'ABCDEFG'
        
        DBSession.add(MS_peptide)
        DBSession.flush()
        
        series_header = ['data:time(min):0','data:time(min):5','data:time(min):20', 'data:stddev(min):5', 'data:stddev(min):20']
        series = [0,4,1,3,2]
        
        t = tasks.insert_run_data.apply_async((MS_peptide, series_header, "run1", series))
        
        t.get()
        assert t.successful()
        
        DBSession.flush()
        
        result = DBSession.query(experiment.ExperimentData).filter_by(MS_id=MS_peptide.id).all()
        result = sorted(result, key=lambda item: item.priority)
        
        self.assertEqual(5, len(result))
        
        self.assertEqual("run1", result[0].run)
        self.assertEqual("time(min)", result[0].type)
        self.assertEqual('0', result[0].label)
        self.assertEqual(0, result[0].value)
        
        self.assertEqual("run1", result[1].run)
        self.assertEqual("time(min)", result[1].type)
        self.assertEqual('5', result[1].label)
        self.assertEqual(4, result[1].value)

        self.assertEqual("run1", result[2].run)
        self.assertEqual("time(min)", result[2].type)
        self.assertEqual('20', result[2].label)
        self.assertEqual(1, result[2].value)

        self.assertEqual("run1", result[3].run)
        self.assertEqual("stddev(min)", result[3].type)
        self.assertEqual('5', result[3].label)
        self.assertEqual(3, result[3].value)

        self.assertEqual("run1", result[4].run)
        self.assertEqual("stddev(min)", result[4].type)
        self.assertEqual('20', result[4].label)
        self.assertEqual(2, result[4].value)
        
    
    @patch('ptmworker.tasks.finalize_import')
    @patch('ptmworker.tasks.insert_run_data')
    @patch('ptmworker.tasks.load_peptide')
    @patch('ptmworker.tasks.load_proteins')
    def test_start_import_should_generate_subtasks_for_input_file(self, patch_loadProtein, patch_loadPeptide, patch_load_run_data, patch_finalize):
        os.chdir(settings.ptmscout_path)
        
        exp = createMockExperiment()
        exp.datafile = os.path.join("test", "test_dataset_formatted.txt")
        
        col_map = {'accession':0, 'peptide':3, 'run':4, 'data':range(5, 21)}
        res = tasks.start_import.apply_async((exp, col_map))
        
        process_id = res.get()
        assert res.successful()

        assert patch_loadProtein.s.called
        args = ""
        for i in xrange(0, len(patch_loadProtein.s.call_args_list)):
            args += str(patch_loadProtein.s.call_args_list[i])
        
        self.assertEqual(3, len( patch_loadProtein.s.call_args_list ))
        self.assertEqual( args.find('P07197'), args.rfind('P07197')) 
        
        assert call(exp.id, 'P50914', 'AALLKAsPK') in patch_loadPeptide.s.call_args_list
        assert call(exp.id, 'Q8N9T8', 'AFVEDsEDEDGAGEGGSSLLQK') in patch_loadPeptide.s.call_args_list
        assert call(exp.id, 'Q6KC79', 'AITSLLGGGsPK') in patch_loadPeptide.s.call_args_list
        
        self.assertEquals(17, len(patch_load_run_data.s.call_args_list))
        
        patch_finalize.s.assert_called_once_with(exp)
        
        self.assertEqual(process_id, exp.import_process_id)
        self.assertEqual('loading', exp.status)
        exp.saveExperiment.assert_called_once_with()
        
    @patch('ptmworker.tasks.insert_run_data')
    @patch('ptmworker.tasks.load_peptide')
    @patch('ptmworker.tasks.load_proteins')
    def test_start_import_should_generate_subtasks_for_input_file_should_work_when_exceeding_batch_size(self, patch_loadProtein, patch_loadPeptide, patch_load_run_data):
        os.chdir(settings.ptmscout_path)
        
        exp = createMockExperiment()
        exp.datafile = os.path.join("test", "test_dataset_formatted.txt")
        
        col_map = {'accession':0, 'peptide':3, 'run':4, 'data':range(5, 21)}
        res = tasks.start_import.apply_async((exp, col_map, 5))
        
        process_id = res.get()
        assert res.successful()
        
        self.assertEqual(3, len( patch_loadProtein.s.call_args_list ))