from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.database import upload

class TestConfigureAnnotationUploadView(UnitTestCase):
    def test_view(self):
        pass
        
class IntegrationTestConfigureAnnotationUploadView(IntegrationTestCase):
    def test_view_integration(self):
        session = upload.Session()
        session.user_id = self.bot.user.id
        session.change_description = ""
        session.data_file = 'test/experiment.28.test.annotations.tsv'
        session.load_type = 'annotations'
        session.stage='config'
        session.save()
        
        result = self.ptmscoutapp.get("/experiments/28/annotate/%d/configure" % (session.id), status=200)
        result.form.submit(status=302)