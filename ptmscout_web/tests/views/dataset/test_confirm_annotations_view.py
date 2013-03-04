from tests.PTMScoutTestCase import IntegrationTestCase
from ptmscout.config import strings
from mock import patch
from ptmscout.database import upload
        
class IntegrationTestUploadStatusView(IntegrationTestCase):
    def test_view_upload_already_started(self):
        self.bot.login()
        
        session = upload.Session()
        session.experiment_id = 1
        session.user_id = self.bot.user.id
        session.load_type='new'
        session.data_file='exp_file'
        session.parent_experiment=None
        session.change_description=''
        session.stage='complete'
        session.save()
        
        result = self.ptmscoutapp.get("/experiments/28/annotate/%d/confirm" % (session.id), status=200)
        result.mustcontain(strings.annotation_upload_started_page_title)
    
    @patch('ptmworker.annotate_tasks.start_annotation_import.apply_async')
    def test_view_integration(self, patch_startImport):
        from ptmscout.database import experiment
        self.bot.login()
        
        exp = experiment.getExperimentById(1, 0, False)
        exp.status = 'configuration'
        exp.saveExperiment()
        
        session = upload.Session()
        session.experiment_id = 1
        session.user_id = self.bot.user.id
        session.load_type='new'
        session.data_file='exp_file'
        session.parent_experiment=None
        session.change_description=''
        session.stage='confirm'
        session.save()
        
        result = self.ptmscoutapp.get("/experiments/28/annotate/%d/confirm" % session.id, status=200)
                
        result.mustcontain(strings.annotation_upload_confirm_message)
        result.form.set('terms_of_use', True)
        result = result.form.submit()
        
        result.mustcontain(strings.annotation_upload_started_message % ("http://localhost/account/experiments"))
        
        assert patch_startImport.called
