from mock import patch
from ptmscout.config import strings, settings
from tests.PTMScoutTestCase import IntegrationTestCase
import os
        
class IntegrationTestMCAMEnrichmentView(IntegrationTestCase):
    @patch('ptmworker.mcam_tasks.run_mcam_analysis.apply_async')
    def test_view_integration(self, patch_run_mcam):
        result = self.ptmscoutapp.get("/experiments/26/mcam_enrichment", status=200)
        result.mustcontain(strings.mcam_enrichment_page_title)
        
        d = result.pyquery
        self.assertIn("Number of Clusters: 0", d('.grid_3').text())
        
        result = result.form.submit().follow()
        result.mustcontain(strings.mcam_enrichment_started_page_title)
        
        self.assertTrue( patch_run_mcam.called )
        
        
    def test_view_download_should_function(self):
        user_id = self.bot.user.id
        mcam_id = 1234567890
        exp_id = 26
        
        mcam_filename = "%d.mcam.%d.%d" % (exp_id, user_id, mcam_id)
        mcam_path = os.path.join(settings.ptmscout_path, settings.mcam_file_path, mcam_filename)
        f = open(mcam_path, 'w')
        f.write("These are the file contents")
        f.close()
        
        result = self.ptmscoutapp.get("/experiments/%d/mcam_download/%d" % (exp_id, mcam_id))
        
        result.mustcontain("These are the file contents")
        
        os.remove(mcam_path)