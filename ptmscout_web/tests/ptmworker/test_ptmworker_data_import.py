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
    @patch('ptmworker.tasks.load_protein')
    def test_start_import_should_generate_subtasks_for_input_file(self, patch_loadProtein, patch_loadPeptide, patch_load_run_data):
        os.chdir(settings.ptmscout_path)
        
        exp = createMockExperiment()
        exp.datafile = os.path.join("test", "test_dataset.txt")
        
        col_map = {'accession':0, 'peptide':2, 'run':3, 'data':range(4, 20)}
        res = tasks.start_import.apply_async((exp, col_map))
        
        process_id = res.get()
        assert res.successful()

        assert call('P07197') in patch_loadProtein.s.call_args_list
        assert call('P27824') in patch_loadProtein.s.call_args_list
        assert call('P50914') in patch_loadProtein.s.call_args_list
        
        assert call('AALLKApSPK') in patch_loadPeptide.s.call_args_list
        assert call('AFVEDpSEDEDGAGEGGSSLLQK') in patch_loadPeptide.s.call_args_list
        assert call('AITSLLGGGpSPK') in patch_loadPeptide.s.call_args_list
        assert call('ELSNSPLRENpSFGSPLEFR') in patch_loadPeptide.s.call_args_list
        
        self.assertEquals(18, len(patch_load_run_data.s.call_args_list))
        
        self.assertEqual(1, len([ c for c in patch_loadProtein.s.call_args_list if c == call('P07197')]))
        
        self.assertEqual(process_id, exp.import_process_id)
        exp.saveExperiment.assert_called_once_with()