from ptmscout.views.experiment.download_view import download_experiment
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockExperiment, createMockUser,\
    createMockError
from mock import patch
from pyramid.httpexceptions import HTTPForbidden
from ptmscout.utils import uploadutils
from ptmscout.config import strings
from ptmscout.database import experiment
import os, re
from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase

class TestExperimentDownloadView(UnitTestCase):
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_download_experiment_should_raise_forbidden_to_non_owner(self, patch_getExperiment):
        data_file = 'test/test_dataset_formatted.txt'
        exp = createMockExperiment()
        exp.dataset = data_file
        
        createMockError(1, "1 Error", experiment=exp)
        createMockError(8, "2 Error", experiment=exp)
        createMockError(8, "3 Error", experiment=exp)
        createMockError(17, "4 Error", experiment=exp)
        
        user = createMockUser()
        user.myExperiments.return_value = []
        
        request = DummyRequest()
        request.matchdict['id'] = exp.id
        request.user = user
        
        patch_getExperiment.return_value = exp
        
        try:
            download_experiment(request)
        except HTTPForbidden:
            pass
        else:
            self.fail("Expected HTTPForbidden")
        
        patch_getExperiment.assert_called_once_with(exp.id, user)
        user.myExperiments.assert_called_once_with()
    
        
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_download_experiment_view_result_should_contain_recorded_errors(self, patch_getExperiment):
        data_file = 'test/test_dataset_formatted.txt'
        exp = createMockExperiment()
        exp.dataset = data_file
        
        exp_headers, exp_data = uploadutils.load_header_and_data_rows(data_file, 100)
        exp_headers.insert(0, strings.experiment_upload_error_reasons_column_title)
        
        for i in xrange(0, len(exp_data)):
            exp_data[i].insert(0, "")
            
        exp_data[0][0] = "1 Error"
        exp_data[7][0] = "2 Error, 3 Error"
        exp_data[16][0] = "4 Error"
        
        exp_data = [row for row in exp_data if row[0] != ""]
        
        createMockError(1, "1 Error", experiment=exp)
        createMockError(8, "2 Error", experiment=exp)
        createMockError(8, "3 Error", experiment=exp)
        createMockError(17, "4 Error", experiment=exp)
        
        user = createMockUser()
        user.myExperiments.return_value = [exp]
        
        request = DummyRequest()
        request.GET['errors'] = "True"
        request.matchdict['id'] = exp.id
        request.user = user
        
        patch_getExperiment.return_value = exp
        
        result = download_experiment(request)
        
        patch_getExperiment.assert_called_once_with(exp.id, user)
        user.myExperiments.assert_called_once_with()
        
        
        self.assertEqual(exp_headers, result['header'])
        self.assertEqual(exp_data, result['data'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_download_experiment_view(self, patch_getExperiment):
        data_file = 'test/test_dataset_formatted.txt'
        exp = createMockExperiment(eid=20)
        exp.dataset = data_file
        
        exp_headers, exp_data = uploadutils.load_header_and_data_rows(data_file, 100)
        
        createMockError(1, "1 Error", experiment=exp)
        createMockError(8, "2 Error", experiment=exp)
        createMockError(8, "3 Error", experiment=exp)
        createMockError(17, "4 Error", experiment=exp)
        
        user = createMockUser()
        user.myExperiments.return_value = [exp]
        
        request = DummyRequest()
        request.matchdict['id'] = exp.id
        request.user = user
        
        patch_getExperiment.return_value = exp
        
        result = download_experiment(request)
        
        patch_getExperiment.assert_called_once_with(exp.id, user)
        user.myExperiments.assert_called_once_with()
        
        self.assertEqual(exp_headers, result['header'])
        self.assertEqual(exp_data, result['data'])
        self.assertEqual('text/tab-separated-values', request.response.content_type)
        self.assertEqual('attachment; filename="experiment.20.tsv"', request.response.content_disposition)
        
        
class IntegrationTestExperimentDownloadView(IntegrationTestCase):
    def test_view_should_forbidden_if_not_logged_in(self):
        self.bot.logout()
        result = self.ptmscoutapp.get("/experiments/28/download", status=403)
        result.mustcontain("forbidden")
        
    def test_view_should_forbidden_if_not_owner(self):
        self.bot.login()
        result = self.ptmscoutapp.get("/experiments/28/download", status=403)
        result.mustcontain("forbidden")
        
    def test_view_should_succeed_if_owner(self):
        self.bot.login()
        
        exp = experiment.getExperimentById(28, secure=False)
        exp.dataset = os.path.join('test','test_dataset_formatted.txt')
        exp.saveExperiment()
        
        self.bot.acquire_experiments([28])
        
        self.ptmscoutapp.get("/experiments/28/download", status=200)
        
    def test_export_view_is_private(self):
        self.bot.logout()
        result = self.ptmscoutapp.get("/experiments/28/export", status=403)
        result.mustcontain("forbidden")

    def log_abort(self):
        import logging 
        logging.getLogger('ptmscout').debug("Transaction aborted")
        raise Exception()

    def session_flush(self):
        import logging
        from ptmscout.database import DBSession
        DBSession.flush()
        logging.getLogger('ptmscout').debug("Transaction committed")

    @patch('transaction.abort')
    @patch('transaction.commit')
    @patch('ptmscout.utils.mail.celery_send_mail')
    def test_export_view_should_have_default_columns(self, patch_mail, patch_commit, patch_abort):
        self.bot.login()
        patch_commit.side_effect = self.session_flush
        patch_abort.side_effect = self.log_abort

        result = self.ptmscoutapp.get("/experiments/28/export?annotate=no")
        result.mustcontain(strings.experiment_export_started_page_title)

        m = re.search(r'<a href="(.*)">here</a>', str(patch_mail.call_args_list))

        download_url = m.group(1)
        last_path_seg = download_url[download_url.rfind('/')+1:]
        result = self.ptmscoutapp.get(download_url)
                
        exp_header = 'MS_id    query_accession    acc_gene	locus protein_name    species    peptide    mod_sites    gene_site    aligned_peptides    modification_types    24H_EGF:data:time:0    24H_EGF:data:time:5    24H_EGF:data:time:10    24H_EGF:data:time:30    24H_EGF:stddev:time:0    24H_EGF:stddev:time:10    24H_EGF:stddev:time:30    24H_HRG:data:time:0    24H_HRG:data:time:5    24H_HRG:data:time:10    24H_HRG:data:time:30    24H_HRG:stddev:time:0    24H_HRG:stddev:time:10    24H_HRG:stddev:time:30    P_EGF:data:time:0    P_EGF:data:time:5    P_EGF:data:time:10    P_EGF:data:time:30    P_EGF:stddev:time:0    P_EGF:stddev:time:10    P_EGF:stddev:time:30    P_HRG:data:time:0    P_HRG:data:time:5    P_HRG:data:time:10    P_HRG:data:time:30    P_HRG:stddev:time:0    P_HRG:stddev:time:10    P_HRG:stddev:time:30'.replace('    ', '\t')
        exp_header = exp_header.split()
        
        lines = str(result).split("\n")
        self.assertEqual('Response: 200 OK', lines[0])
        self.assertEqual('Content-Disposition: attachment; filename="experiment.28.%d.%s.tsv"' % (self.bot.user.id, last_path_seg), lines[1])
        self.assertEqual('Content-Type: text/tab-separated-values; charset=UTF-8', lines[2])
        
        for line in lines[3:]:
            row = line.split("\t")
            self.assertEqual(len(exp_header), len(row))
 
    @patch('transaction.abort')
    @patch('transaction.commit')
    @patch('ptmscout.utils.mail.celery_send_mail')
    def test_export_view_should_annotate_experiment(self, patch_mail, patch_commit, patch_abort):
        self.bot.login()
        patch_commit.side_effect = self.session_flush
        patch_abort.side_effect = self.log_abort

        result = self.ptmscoutapp.get("/experiments/28/export?annotate=yes")
        result.mustcontain(strings.experiment_export_started_page_title)

        m = re.search(r'<a href="(.*)">here</a>', str(patch_mail.call_args_list))
        
        download_url = m.group(1)
        last_path_seg = download_url[download_url.rfind('/')+1:]
        result = self.ptmscoutapp.get(download_url)
         
        exp_header = 'MS_id    query_accession 	acc_gene	locus protein_name    species    peptide    mod_sites    gene_site    aligned_peptides    modification_types    24H_EGF:data:time:0    24H_EGF:data:time:5    24H_EGF:data:time:10    24H_EGF:data:time:30    24H_EGF:stddev:time:0    24H_EGF:stddev:time:10    24H_EGF:stddev:time:30    24H_HRG:data:time:0    24H_HRG:data:time:5    24H_HRG:data:time:10    24H_HRG:data:time:30    24H_HRG:stddev:time:0    24H_HRG:stddev:time:10    24H_HRG:stddev:time:30    P_EGF:data:time:0    P_EGF:data:time:5    P_EGF:data:time:10    P_EGF:data:time:30    P_EGF:stddev:time:0    P_EGF:stddev:time:10    P_EGF:stddev:time:30    P_HRG:data:time:0    P_HRG:data:time:5    P_HRG:data:time:10    P_HRG:data:time:30    P_HRG:stddev:time:0    P_HRG:stddev:time:10    P_HRG:stddev:time:30    scansite_bind    scansite_kinase    nearby_modifications    nearby_mutations    site_domains    site_regions    protein_domains    protein_GO_BP    protein_GO_CC    protein_GO_MF'.replace('    ', '\t')
        exp_header = exp_header.split()
        
        lines = str(result).split("\n")
        self.assertEqual('Response: 200 OK', lines[0])
        self.assertEqual('Content-Disposition: attachment; filename="experiment.28.%d.%s.tsv"' % (self.bot.user.id, last_path_seg), lines[1])
        self.assertEqual('Content-Type: text/tab-separated-values; charset=UTF-8', lines[2])
        
        for line in lines[3:]:
            row = line.split("\t")
            self.assertEqual(len(exp_header), len(row))
        
