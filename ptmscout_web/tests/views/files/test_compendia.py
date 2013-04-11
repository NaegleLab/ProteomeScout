from pyramid.testing import DummyRequest
from ptmscout.config import strings
from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.views.files import compendia
from mock import Mock

class TestCompendiaDownloadIntegration(IntegrationTestCase):
    def test_compendia_list_should_fail_for_non_users(self):
        self.bot.logout()
        result = self.ptmscoutapp.get('/compendia', status=403)
        result.mustcontain('Forbidden')

    def test_compendia_download_should_fail_for_non_users(self):
        self.bot.logout()
        result = self.ptmscoutapp.get('/compendia/vertebrata.tsv', status=403)
        result.mustcontain('Forbidden')

    def test_compendia_list_should_show_available_compendia(self):
        result = self.ptmscoutapp.get('/compendia', status=200)
        result.mustcontain('everything.tsv')
        result.mustcontain('All proteins and modifications')
        result.mustcontain('ubiquitination.tsv')
        result.mustcontain('glycosylation.tsv')

    def test_compendia_download_should_raise_404(self):
        self.ptmscoutapp.get('/compendia/not_a_real_file.tsv', status=404)

    def test_compendia_download_should_get_file(self):
        result = self.ptmscoutapp.get('/compendia/vertebrata.tsv', status=200)
        result.mustcontain('accessions	acc_gene	locus	protein_name	species	sequence	modifications	domains	mutations	scansite_predictions	GO_terms')

class TestCompendiaDownload(UnitTestCase):
    def test_compendia_list(self):
        request = DummyRequest()
        request.route_url = Mock()
        request.route_url.return_value = 'link'
        result = compendia.compendia_listing(request)

        self.assertEqual(strings.compendia_download_page_title, result['pageTitle'])
        self.assertEqual(strings.compendia_download_page_desc, result['desc'])
        self.assertEqual(8, len(result['files']))
