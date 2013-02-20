from ptmworker import notify_tasks
from tests.PTMScoutTestCase import IntegrationTestCase
from tests.views.mocking import createMockExperiment
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
    @patch('ptmscout.database.experiment.setExperimentProgress')
    def test_set_progress(self, patch_setProgress, patch_commit):
        notify_tasks.set_progress(2000, 100, 1000)
        patch_setProgress.assert_called_once_with(2000, 100, 1000)
        patch_commit.assert_called_once_with()


    @patch('transaction.commit')
    @patch('ptmscout.database.experiment.getExperimentById')
    @patch('ptmscout.database.experiment.setExperimentProgress')
    def test_set_loading_stage(self, patch_setProgress, patch_getExp, patch_commit):
        exp = createMockExperiment()
        exp.loading_stage = 'query'
        patch_getExp.return_value = exp
        
        notify_tasks.set_loading_stage(exp.id, 'proteins', 1000)
        patch_getExp.assert_called_once_with(exp.id, check_ready=False, secure=False)
        
        self.assertEqual('proteins', exp.loading_stage)
        exp.saveExperiment.assert_called_once_with()

        patch_setProgress.assert_called_once_with(exp.id, 0, 1000)

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
        notify_tasks.finalize_experiment_error_state_callback(uuid, exp_id, user_email, application_url)

        patch_result.assert_called_once_with(uuid)

        patch_getExp.assert_called_once_with(exp.id, check_ready=False, secure=False)
        self.assertEqual('some traceback', exp.failure_reason)
        self.assertEqual('error', exp.status)
        exp.saveExperiment.assert_called_once_with()

        email_subject = strings.experiment_upload_failed_subject
        email_message = strings.experiment_upload_failed_message % (exp.name, exp.loading_stage, "Exception: Disturbance in the force", application_url)
        patch_sendmail.assert_called_once_with([user_email, settings.adminEmail], email_subject, email_message)

        patch_commit.called


    @patch('ptmscout.utils.mail.celery_send_mail')
    @patch('transaction.commit')
    @patch('ptmscout.database.experiment.countErrorsForExperiment')
    @patch('ptmscout.database.modifications.countMeasuredPeptidesForExperiment')
    @patch('ptmscout.database.modifications.countProteinsForExperiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_finalize_import(self, patch_getExp, patch_countProteins, patch_countPeptides, patch_countErrors, patch_commit, patch_sendmail):
        exp = createMockExperiment()
        exp.status = 'loading'

        pep_cnt = 10
        prot_cnt = 7
        err_cnt = 1
        patch_countProteins.return_value=prot_cnt
        patch_countPeptides.return_value=pep_cnt
        patch_countErrors.return_value=err_cnt

        exp_id = exp.id
        user_email = "user@institute.edu"
        application_url = "http://ptmscout.example.com"

        patch_getExp.return_value = exp
        notify_tasks.finalize_import(exp_id, user_email, application_url)

        patch_countErrors.assert_called_once_with(exp_id)
        patch_countPeptides.assert_called_once_with(exp_id)
        patch_countProteins.assert_called_once_with(exp_id)

        patch_getExp.assert_called_once_with(exp.id, check_ready=False, secure=False)
        self.assertEqual('loaded', exp.status)
        exp.saveExperiment.assert_called_once_with()

        email_subject = strings.experiment_upload_finished_subject
        error_log_url = "%s/experiments/%d/errors" % (application_url, exp_id)
        email_message = strings.experiment_upload_finished_message % (exp.name, pep_cnt, prot_cnt, err_cnt, error_log_url)
        patch_sendmail.assert_called_once_with([user_email], email_subject, email_message)

        patch_commit.assert_called_once_with()


