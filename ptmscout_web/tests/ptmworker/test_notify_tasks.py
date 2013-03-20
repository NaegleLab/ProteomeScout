from ptmworker import notify_tasks
from tests.PTMScoutTestCase import IntegrationTestCase
from tests.views.mocking import createMockExperiment, createMockJob,\
    createMockUser
from mock import patch
from ptmscout.config import strings

class NotifyTasksTest(IntegrationTestCase):
    @patch('transaction.commit')
    @patch('ptmscout.database.jobs.getJobById')
    def test_set_loading_status_should_get_job_and_mark(self, patch_getJob, patch_commit):
        job = createMockJob()
        job.status = 'loading'
        patch_getJob.return_value = job
        
        notify_tasks.set_job_status(job.id, 'loaded')
        patch_getJob.assert_called_once_with(job.id)
        
        self.assertEqual('loaded', job.status)
        job.save.assert_called_once_with()

        patch_commit.assert_called_once_with()

    @patch('transaction.commit')
    @patch('ptmscout.database.jobs.getJobById')
    def test_set_progress(self, patch_getJob, patch_commit):
        job = createMockJob()
        job.status = 'loading'
        job.progress = 0
        job.max_progess = 0
        patch_getJob.return_value = job
        
        notify_tasks.set_job_progress(2000, 100, 1000)
        patch_getJob.assert_called_once_with(2000)
        
        self.assertEqual(100, job.progress)
        self.assertEqual(1000, job.max_progress)
        
        patch_commit.assert_called_once_with()


    @patch('transaction.commit')
    @patch('ptmscout.database.jobs.getJobById')
    def test_set_loading_stage(self, patch_getJob, patch_commit):
        job = createMockJob()
        job.status = 'loading'
        job.stage = 'query'
        job.progress = 0
        job.max_progess = 0
        
        patch_getJob.return_value = job
        
        notify_tasks.set_job_stage(job.id, 'proteins', 1000)
        patch_getJob.assert_called_once_with(job.id)
        
        self.assertEqual('proteins', job.stage)
        self.assertEqual(0, job.progress)
        self.assertEqual(1000, job.max_progress)
        job.save.assert_called_once_with()

        patch_commit.assert_called_once_with()


    @patch('ptmscout.utils.mail.celery_send_mail')
    @patch('transaction.commit')
    @patch('ptmscout.database.experiment.countErrorsForExperiment')
    @patch('ptmscout.database.modifications.countMeasuredPeptidesForExperiment')
    @patch('ptmscout.database.modifications.countProteinsForExperiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_finalize_import(self, patch_getExp, patch_countProteins, patch_countPeptides, patch_countErrors, patch_commit, patch_sendmail):
        exp = createMockExperiment()
        exp.status = 'loading'
        exp.job = createMockJob()
        
        application_url = "http://ptmscout.example.com"
        user_email = "user@institute.edu"
        exp_id = exp.id
        
        exp.job.user = createMockUser(email=user_email)
        
        exp.job.result_url = "%s/experiments/%d" % (application_url, exp_id)

        pep_cnt = 10
        prot_cnt = 7
        err_cnt = 1
        patch_countProteins.return_value=prot_cnt
        patch_countPeptides.return_value=pep_cnt
        patch_countErrors.return_value=err_cnt


        patch_getExp.return_value = exp
        notify_tasks.finalize_experiment_import(exp_id)

        patch_countErrors.assert_called_once_with(exp_id)
        patch_countPeptides.assert_called_once_with(exp_id)
        patch_countProteins.assert_called_once_with(exp_id)

        patch_getExp.assert_called_once_with(exp.id, check_ready=False, secure=False)
        
        exp.job.finish.assert_called_once_with()
        exp.job.save.assert_called_once_with()

        email_subject = strings.experiment_upload_finished_subject
        error_log_url = "%s/errors" % (exp.job.result_url)
        email_message = strings.experiment_upload_finished_message % (exp.name, pep_cnt, prot_cnt, err_cnt, error_log_url)
        patch_sendmail.assert_called_once_with([user_email], email_subject, email_message)

        patch_commit.assert_called_once_with()


