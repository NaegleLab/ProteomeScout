from tests.PTMScoutTestCase import IntegrationTestCase
from mock import patch, call
from ptmworker import tasks
import os
from ptmscout.database import experiment
from ptmscout.config import settings
from tests.views.mocking import createMockExperiment

class PTMWorkDataImportTestCase(IntegrationTestCase):
    
    @patch('ptmworker.tasks.insert_run_data')
    @patch('ptmworker.tasks.load_peptide')
    @patch('ptmworker.tasks.load_proteins')
    def test_start_import_should_generate_subtasks_for_input_file(self, patch_loadProtein, patch_loadPeptide, patch_load_run_data):
        os.chdir(settings.ptmscout_path)
        
        exp = createMockExperiment()
        exp.datafile = os.path.join("test", "test_dataset.txt")
        
        col_map = {'accession':0, 'peptide':2, 'run':3, 'data':range(4, 20)}
        res = tasks.start_import.apply_async((exp, col_map))
        
        process_id = res.get()
        assert res.successful()

        assert patch_loadProtein.s.called
        args = ""
        for i in xrange(0, len(patch_loadProtein.s.call_args_list)):
            args += str(patch_loadProtein.s.call_args_list[i])
        self.assertEqual( args.find('P07197'), args.rfind('P07197')) 
        
        assert call('P50914', 'AALLKApSPK') in patch_loadPeptide.s.call_args_list
        assert call('Q8N9T8', 'AFVEDpSEDEDGAGEGGSSLLQK') in patch_loadPeptide.s.call_args_list
        assert call('Q6KC79', 'AITSLLGGGpSPK') in patch_loadPeptide.s.call_args_list
        assert call('A0AUK8', 'ELSNSPLRENpSFGSPLEFR') in patch_loadPeptide.s.call_args_list
        
        self.assertEquals(18, len(patch_load_run_data.s.call_args_list))
        
        self.assertEqual(process_id, exp.import_process_id)
        self.assertEqual('loading', exp.status)
        exp.saveExperiment.assert_called_once_with()