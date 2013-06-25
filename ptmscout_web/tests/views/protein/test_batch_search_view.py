from mock import patch
from ptmscout.config import strings
from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase
import StringIO, re, csv

class TestProteinSearchViewIntegration(IntegrationTestCase):

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
    def test_integration(self, patch_mail, patch_commit, patch_abort):
        patch_commit.side_effect = self.session_flush
        patch_abort.side_effect = self.log_abort

        result = self.ptmscoutapp.get('/batch', status=200)
        result.mustcontain( strings.protein_batch_search_page_title )

        acc_list = ['P00533', 'P04626', 'P53663']

        result.form.set('terms_of_use', True)
        result.form.set('accessions', '\n'.join(acc_list))
        result = result.form.submit().follow()
        result.mustcontain( strings.protein_batch_search_submitted_page_title )

        email = str(patch_mail.call_args_list)
        m = re.search(r'<a href="(.*)">here</a>', email)
        assert m != None, "Download link not found in email: " + email

        assert 'Errors Encountered: 1' in email, email
        assert 'Proteins Loaded: 2' in email, email

        result = self.ptmscoutapp.get(m.group(1), status=200)

        f = StringIO.StringIO(result.body)
        dr = csv.DictReader(f, dialect='excel-tab')

        error_cnt = 0
        for row in dr:
            assert row['query_accession'] in acc_list
            if 'ERRORS' in row['protein_id']:
                error_cnt += 1

        assert error_cnt == 1
