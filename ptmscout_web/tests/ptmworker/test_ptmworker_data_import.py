# Skip
from tests.PTMScoutTestCase import IntegrationTestCase
from mock import patch
from ptmworker import tasks
from ptmscout.config import strings
from tests.views.mocking import createMockExperiment, createMockMeasurement,\
    createMockError, createMockSession, createMockUser, createMockProtein,\
    createMockPTM, createMockPhosphopep
from ptmscout.utils import uploadutils

class PTMWorkDataImportTestCase(IntegrationTestCase):

    @patch('ptmscout.utils.mail.celery_send_mail')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')
    @patch('ptmworker.upload_helpers.mark_experiment')
    def test_finalize_import_should_set_experiment_to_loaded_status_and_statistics(self, patch_mark_experiment, patch_getPeptides, patch_send_mail):
        exp = createMockExperiment()
        e1 = createMockError(2, "Some error", experiment = exp)
        patch_mark_experiment.return_value = exp
        
        m1 = createMockMeasurement(10, exp.id)
        m2 = createMockMeasurement(10, exp.id)
        m3 = createMockMeasurement(12, exp.id)
        
        patch_getPeptides.return_value = [m1,m2,m3]
        
        user_email = 'someguy@institute.edu'
        app_url = 'http://example.com'
        
        res = tasks.finalize_import.apply_async((exp.id, user_email, app_url))
        res.get()
        assert res.successful()
        
        patch_mark_experiment.assert_called_once_with(exp.id, 'loaded')
        patch_getPeptides.assert_called_once_with(exp.id, check_ready=False, secure=False)
        
        expected_message = strings.experiment_upload_finished_message % (exp.name, 3, 2, 1, app_url + "/experiments/%d/errors" % (exp.id))
        patch_send_mail.assert_called_once_with([user_email], strings.experiment_upload_finished_subject, expected_message)

    @patch('ptmscout.database.experiment.createExperimentError')
    @patch('ptmworker.upload_helpers.create_modifications')
    def test_load_peptide_should_aggregate_errors(self, patch_create_mods, patch_create_error):
        exp_id = 10
        patch_create_mods.side_effect = uploadutils.ParseError(None, None, "An error")
        
        acc = "Q3K102"
        protein_id = 10
        prot_seq = "ABSFGWERJADSFKJ"
        taxonomy = ["Bacteria", "Escherichia"]
        pep_seq = "GWeRJAD"
        modlist = [createMockPTM()]
        run_task_args = [(1, "a","b","c","d"),(3, "d","e","f","g")]
        line_mapping = {1: (acc, pep_seq),
                        3: (acc, pep_seq)}
        
        res = tasks.load_peptide.apply_async(((protein_id, prot_seq, taxonomy), exp_id, pep_seq, modlist, line_mapping, run_task_args))
        res.get()
        
        self.assertTrue(res.successful())
        patch_create_error.assert_any_call(exp_id, 1, acc, pep_seq, "An error")
        patch_create_error.assert_any_call(exp_id, 3, acc, pep_seq, "An error")
        

    @patch('ptmworker.upload_helpers.insert_run_data')
    @patch('ptmscout.database.modifications.MeasuredPeptide')
    @patch('ptmworker.upload_helpers.create_modifications')
    def test_load_peptide_should_measured_peptide_record_with_modifications_and_insert_run_data(self, patch_create_mods, patch_measurement, patch_insert_run):
        exp_id = 10
        prot = createMockProtein()
        pep1 = createMockPhosphopep(prot.id)
        pep2 = createMockPhosphopep(prot.id)
        
        patch_create_mods.return_value = [pep1, pep2]
        measuredPeptide = patch_measurement.return_value
        measuredPeptide.phosphopeps = []
        
        protein_id = prot.id
        pep_seq = "GWeRjAD"
        
        line_mapping = {1:("Q3K102", pep_seq), 3:("Q3K102", pep_seq)}
        prot_seq = "ABSFGWERJADSFKJ"
        taxonomy = ["Bacteria", "Escherichia"]
        modlist = [createMockPTM(),createMockPTM()]
        run_task_args = [(1, "a","b","c","d"),(3, "d","e","f","g")]
        
        res = tasks.load_peptide.apply_async(((protein_id, prot_seq, taxonomy), exp_id, pep_seq, modlist, line_mapping, run_task_args))
        res.get()
        
        self.assertTrue(res.successful())
        
        self.assertEqual(exp_id, measuredPeptide.experiment_id)
        self.assertEqual(pep_seq, measuredPeptide.phosphopep)
        self.assertEqual(prot.id, measuredPeptide.protein_id)
        self.assertEqual([pep1,pep2], measuredPeptide.phosphopeps)
        
        measuredPeptide.save.assert_called_once_with()
        
        patch_insert_run.assert_any_call(measuredPeptide, 1, "a","b","c","d")
        patch_insert_run.assert_any_call(measuredPeptide, 3, "d","e","f","g")
        
        

    @patch('ptmscout.database.experiment.createExperimentError')
    @patch('ptmworker.upload_helpers.find_protein')
    def test_load_protein_should_aggregate_errors(self, patch_find_prot, patch_create_error):
        exp_id = 10
        prot = createMockProtein()
        patch_find_prot.side_effect = uploadutils.ParseError(None, None, "An error")
        taxons = ["some taxons"]
        
        accessions = ["some", "Accessions"]
        affected_lines = [1,3]
        line_mapping = {1: ('acc', 'pep1'), 3: ('acc', 'pep2')}
        res = tasks.load_protein.apply_async(((prot.name, prot.acc_gene, taxons, prot.species.name, accessions, prot.sequence), exp_id, affected_lines, line_mapping))
        ret_val = res.get()
        
        self.assertTrue(res.successful())
        self.assertEqual(None, ret_val)
        patch_create_error.assert_any_call(exp_id, 1, 'acc', 'pep1', "An error")
        patch_create_error.assert_any_call(exp_id, 3, 'acc', 'pep2', "An error")

    @patch('ptmworker.upload_helpers.find_protein')
    def test_load_protein_should_return_new_protein_info(self, patch_find_prot):
        exp_id = 10
        prot = createMockProtein()
        patch_find_prot.return_value = prot
        taxons = ["some taxons"]
        
        accessions = ["some", "Accessions"]
        line_mapping = "some line map"
        
        res = tasks.load_protein.apply_async(((prot.name, prot.acc_gene, taxons, prot.species.name, accessions, prot.sequence), exp_id, [1,3], line_mapping))
        ret_val = res.get()
        
        self.assertEqual((prot.id, prot.sequence, taxons), ret_val)
        patch_find_prot.assert_called_once_with(prot.name, prot.acc_gene, prot.sequence, accessions, prot.species.name)

    @patch('ptmworker.upload_helpers.mark_experiment')
    def test_process_error_state_should_change_experiment_status(self, patch_mark):
        exp = createMockExperiment()
        patch_mark.return_value = exp
        
        res = tasks.process_error_state.apply_async((exp.id,))
        res.get()
        
        patch_mark.assert_called_once_with(exp.id, 'error')
    
    @patch('ptmworker.tasks.load_peptide')
    @patch('ptmworker.tasks.load_protein')
    def test_create_import_tasks(self, patch_protein, patch_peptide):
        exp_id = 10
        prot_map = {"Q06FX4":"prot info 1", "A0FGVD":"prot info 2"}
        accessions = {"Q06FX4":[1,2,4,5], "A0FGVD":[3,6]}
        line_mapping = {1:("Q06FX4","ABD"),
                        2:("Q06FX4","DEF"),
                        3:("A0FGVD","GHI"),
                        4:("Q06FX4","ABD"),
                        5:("Q06FX4","DEF"),
                        6:("A0FGVD","GHI")
                        }
        
        peptides = {"Q06FX4":["ABD", "DEF"], "A0FGVD":["GHI"]}
        mod_map = {("Q06FX4", "ABD"): "phos", ("Q06FX4", "DEF"): "methylation", ("A0FGVD", "GHI"): "acetylation"}
        series_headers = ["some", "headers"]
        units = "time(min)"
        data_runs = {("Q06FX4", "ABD"): {"run1":(1, [1,2]), "run2":(4, [2,3])}, 
                     ("Q06FX4", "DEF"): {"run1":(2, [3,4]), "run2":(5, [7,8])}, 
                     ("A0FGVD", "GHI"): {"run1":(3, [5,6]), "run2":(6, [9,10])}}
        
        
        import_tasks = tasks.create_import_tasks(exp_id, prot_map, accessions, peptides, mod_map, line_mapping, series_headers, units, data_runs)
        
        patch_peptide.s.assert_any_call(exp_id, "ABD", "phos", line_mapping, [(1, 'time(min)', ["some", "headers"], "run1", [1,2]),(4, 'time(min)', ["some", "headers"], "run2", [2,3])])
        patch_peptide.s.assert_any_call(exp_id, "DEF", "methylation", line_mapping, [(2, 'time(min)', ["some", "headers"], "run1", [3,4]),(5, 'time(min)', ["some", "headers"], "run2", [7,8])])
        patch_peptide.s.assert_any_call(exp_id, "GHI", "acetylation", line_mapping, [(3, 'time(min)', ["some", "headers"], "run1", [5,6]),(6, 'time(min)', ["some", "headers"], "run2", [9,10])])
        
        patch_protein.s.assert_any_call("prot info 1", exp_id, [1,2,4,5], line_mapping)
        patch_protein.s.assert_any_call("prot info 2", exp_id, [3,6], line_mapping)
        
        self.assertEqual(2, len(import_tasks))
        
    
    @patch('ptmscout.database.experiment.createExperimentError')
    @patch('ptmscout.database.upload.getSessionById')
    @patch('ptmworker.tasks.invoke')
    @patch('ptmworker.tasks.finalize_import')
    @patch('ptmworker.tasks.create_import_tasks')
    @patch('ptmworker.upload_helpers.get_proteins_from_ncbi')
    @patch('ptmworker.upload_helpers.get_series_headers')
    @patch('ptmworker.upload_helpers.parse_datafile')
    def test_start_import_should_load_file_get_protein_records_log_errors_and_invoke_subtasks(self, patch_parse, patch_get_headers, patch_get_proteins, patch_create_tasks, patch_finalize, patch_invoke, patch_getSession, patch_createError):
        e1 = uploadutils.ParseError(1, None, "an error")
        e2 = uploadutils.ParseError(9, None, "another error")
        session_id = 100
        exp_id = 2
        app_url = 'http://example.com'
        user = createMockUser()
        session = createMockSession(user, experiment_id = exp_id)
        
        line_mapping = {1: ('Q1G345', "PEP1"),
                        5: ('Q1G345', "PEP2"),
                        10: ('A34PDF', "PEP4"),
                        9: ('A34PDF', "PEP3")
                        }
        accessions = {'Q1G345':[1,5], 'A34PDF':[10,9]}
        patch_parse.return_value = accessions, "some peps", "some mods", "some data", [e1], line_mapping
        patch_get_headers.return_value = "some headers"
        prot_map = {"some accessions":"some records"}
        patch_get_proteins.return_value = prot_map, [e2]
        
        patch_getSession.return_value = session
        
        patch_create_tasks.return_value = "some tasks"
        
        res = tasks.start_import.apply_async((exp_id, session_id, user.email, app_url))
        res.get()
        assert res.successful()
        
        patch_parse.assert_called_once_with(session)
        patch_get_headers.assert_called_once_with(session)
        
        patch_get_proteins.assert_called_once_with(accessions, 1000)
        patch_createError.assert_any_call(exp_id, 1, 'Q1G345', "PEP1", "an error")
        patch_createError.assert_any_call(exp_id, 9, 'A34PDF', "PEP3", "another error")
        
        patch_create_tasks.assert_called_once_with(exp_id, prot_map, accessions, "some peps", "some mods", line_mapping, "some headers", session.units, "some data")
        patch_invoke.assert_called_once_with("some tasks", exp_id, user.email, app_url)
        