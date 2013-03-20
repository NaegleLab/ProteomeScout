from ptmworker import peptide_tasks
from tests.PTMScoutTestCase import IntegrationTestCase
from mock import patch, call

class PeptideTasksTest(IntegrationTestCase):
    @patch('transaction.commit')
    @patch('ptmworker.helpers.upload_helpers.store_stage_input')
    @patch('ptmworker.notify_tasks.set_job_progress')
    @patch('ptmworker.notify_tasks.set_job_stage')
    @patch('ptmworker.peptide_tasks.load_peptide_modification')
    def test_run_peptide_import_should_invoke_load_mod(self, patch_load_peptide, patch_set_stage, patch_set_progress, patch_store_input, patch_commit):
        def save_db():
            from ptmscout.database import DBSession
            DBSession.flush()
        
        prot_map = {"Q9NQG7": 9, "Q9WTS6": 11, "P49848":24, "O88384":105}
        peptides = {"P00533": ["GENETICS"], "Q9NQG7": ["AABDEF", "MEEVFG"], "Q9WTS6": ["GGGTTTCCC"], "P49848": ["LKLKLKLKLK"]}
        exp_id=2600
        job_id=20
        mod_map = { ("P00533", "GENETICS"): ["METHYLATION"],
                    ("Q9NQG7", "AABDEF"): ["DIMETHYLATION"],
                    ("Q9NQG7", "MEEVFG"): ["PHOS"],
                    ("Q9WTS6", "GGGTTTCCC"):["PHOSPHORYLATION"],
                    ("P49848", "LKLKLKLKLK"): ["O-GlcNAc"] }

        data_runs = {("P00533", "GENETICS", "METHYLATION"): {'average': (1, [1,2,3])},
                ("Q9NQG7", "AABDEF", "DIMETHYLATION"): {'average': (1, [1,2,3])},
                ("Q9NQG7", "MEEVFG", "PHOS"): {'average': (1, [1,2,3])},
                ("Q9WTS6", "GGGTTTCCC", "PHOSPHORYLATION"): {'average': (1, [1,2,3])},
                ("P49848", "LKLKLKLKLK", "O-GlcNAc"): {'average': (1, [1,2,3])}}

        patch_commit.side_effect = save_db
        headers = ["some", "sequence", "of", "headers"]
        units = "time(none)"
        load_ambiguities = True


        peptide_tasks.run_peptide_import(prot_map, peptides, mod_map, data_runs, headers, units, load_ambiguities, exp_id, job_id)

        patch_store_input.assert_called_once_with(exp_id, 'peptides', prot_map)
        patch_set_stage.apply_async.assert_called_once_with((job_id, 'peptides', 4))
        patch_commit.assert_called_once_with()
        patch_set_progress.apply_async.assert_called_once_with((job_id, 4, 4))

        self.assertIn(call(2600, True, 'Q9WTS6', 11, 'GGGTTTCCC', 'PHOSPHORYLATION', 'time(none)', ['some', 'sequence', 'of', 'headers'], [(1, 'average', [1, 2, 3])]), patch_load_peptide.call_args_list)
        self.assertIn(call(2600, True, 'Q9NQG7', 9, 'AABDEF', 'DIMETHYLATION', 'time(none)', ['some', 'sequence', 'of', 'headers'], [(1, 'average', [1, 2, 3])]), patch_load_peptide.call_args_list)
        self.assertIn(call(2600, True, 'Q9NQG7', 9, 'MEEVFG', 'PHOS', 'time(none)', ['some', 'sequence', 'of', 'headers'], [(1, 'average', [1, 2, 3])]), patch_load_peptide.call_args_list)
        self.assertIn(call(2600, True, 'P49848', 24, 'LKLKLKLKLK', 'O-GlcNAc', 'time(none)', ['some', 'sequence', 'of', 'headers'], [(1, 'average', [1, 2, 3])]), patch_load_peptide.call_args_list)