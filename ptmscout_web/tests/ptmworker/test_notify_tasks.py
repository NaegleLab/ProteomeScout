from ptmworker import notify_tasks
from tests.PTMScoutTestCase import IntegrationTestCase
from tests.views.mocking import createMockExperiment, createMockError
from mock import patch
from ptmscout.config import settings, strings

class NotifyTasksTest(IntegrationTestCase):
    @patch('transaction.commit')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_set_loading_status_should_get_experiment_and_mark(self, patch_getExp, patch_commit):
        exp = createMockExperiment()
        exp.status = 'loading'
        patch_getExp.return_value = exp
        
        notify_tasks.set_loading_status(exp.id, 'loaded')
        patch_getExp.assert_called_once_with(exp.id, check_ready=False, secure=False)
        
        self.assertEqual('loaded', exp.status)
        exp.saveExperiment.assert_called_once_with()

        patch_commit.assert_called_once_with()

    @patch('transaction.commit')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_set_progress(self, patch_getExp, patch_commit):
        exp = createMockExperiment()
        patch_getExp.return_value = exp
        
        notify_tasks.set_progress(exp.id, 100, 1000)
        patch_getExp.assert_called_once_with(exp.id, check_ready=False, secure=False)
        
        self.assertEqual(100, exp.progress)
        self.assertEqual(1000, exp.max_progress)
        exp.saveExperiment.assert_called_once_with()

        patch_commit.assert_called_once_with()


    @patch('transaction.commit')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_set_loading_stage(self, patch_getExp, patch_commit):
        exp = createMockExperiment()
        patch_getExp.return_value = exp
        
        notify_tasks.set_loading_stage(exp.id, 'proteins', 1000)
        patch_getExp.assert_called_once_with(exp.id, check_ready=False, secure=False)
        
        self.assertEqual('proteins', exp.loading_stage)
        self.assertEqual(0, exp.progress)
        self.assertEqual(1000, exp.max_progress)
        exp.saveExperiment.assert_called_once_with()

        patch_commit.assert_called_once_with()

    @patch('ptmscout.utils.mail.celery_send_mail')
    @patch('celery.result.AsyncResult')
    @patch('transaction.commit')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_finalize_experiment_error_stage(self, patch_getExp, patch_commit, patch_result, patch_sendmail):
        exp = createMockExperiment()
        exp.loading_stage = 'proteins'
        uuid = "some uuid"
        exp_id = exp.id
        user_email = "user@institute.edu"
        application_url = "http://ptmscout.example.com"

        mockResult = patch_result.return_value

        exc = Exception("Disturbance in the force")
        mockResult.traceback = 'some traceback'
        mockResult.get.return_value = exc

        patch_getExp.return_value = exp
        notify_tasks.finalize_experiment_error_state(uuid, exp_id, user_email, application_url)

        patch_result.assert_called_once_with(uuid)

        patch_getExp.assert_called_once_with(exp.id, check_ready=False, secure=False)
        self.assertEqual('some traceback', exp.failure_reason)
        self.assertEqual('error', exp.status)
        exp.saveExperiment.assert_called_once_with()

        email_subject = strings.experiment_upload_failed_subject
        email_message = strings.experiment_upload_failed_message % (exp.name, exp.loading_stage, "Exception: Disturbance in the force", application_url)
        patch_sendmail.assert_called_once_with([user_email, settings.adminEmail], email_subject, email_message)

        patch_commit.assert_called_once_with()


    @patch('ptmscout.utils.mail.celery_send_mail')
    @patch('transaction.commit')
    @patch('ptmscout.database.modifications.countMeasuredPeptidesForExperiment')
    @patch('ptmscout.database.modifications.countProteinsForExperiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_finalize_import(self, patch_getExp, patch_countProteins, patch_countPeptides, patch_commit, patch_sendmail):
        exp = createMockExperiment()
        exp.errors = [createMockError(2, 'failed for peptide')]
        err_cnt = len(exp.errors)
        exp.status = 'loading'

        pep_cnt = 10
        prot_cnt = 7
        patch_countProteins.return_value=prot_cnt
        patch_countPeptides.return_value=pep_cnt

        exp_id = exp.id
        user_email = "user@institute.edu"
        application_url = "http://ptmscout.example.com"

        patch_getExp.return_value = exp
        notify_tasks.finalize_import(exp_id, user_email, application_url)

        patch_getExp.assert_called_once_with(exp.id, check_ready=False, secure=False)
        self.assertEqual('loaded', exp.status)
        exp.saveExperiment.assert_called_once_with()

        email_subject = strings.experiment_upload_finished_subject
        error_log_url = "%s/experiments/%d/errors" % (application_url, exp_id)
        email_message = strings.experiment_upload_finished_message % (exp.name, pep_cnt, prot_cnt, err_cnt, error_log_url)
        patch_sendmail.assert_called_once_with([user_email], email_subject, email_message)

        patch_commit.assert_called_once_with()


