from ptmscout.config import strings
from tests.PTMScoutTestCase import IntegrationTestCase
from ptmscout.database import experiment

class InfoFunctionalTests(IntegrationTestCase):
    def test_about_renderer(self):
        self.ptmscoutapp.get('/about', status=200)
    
    def test_terms_renderer(self):
        self.ptmscoutapp.get('/terms', status=200)
        
    def test_login_views(self):
        self.ptmscoutapp.get('/login')
        self.ptmscoutapp.get('/register')
        
    def test_experiment_not_ready_should_notify_resource_not_ready(self):
        exp = experiment.getExperimentById(1, None)
        exp.job.status = 'running'
        exp.saveExperiment()
        
        response = self.ptmscoutapp.get("/experiments/1")
        response.mustcontain(strings.error_resource_not_ready_page_title)
    
    def test_no_such_session_should_invoke_forbidden(self):
        response = self.ptmscoutapp.get('/upload/2000123/config')
        response.mustcontain("forbidden")
        response = self.ptmscoutapp.get('/upload/2000123/metadata')
        response.mustcontain("forbidden")
        response = self.ptmscoutapp.get('/upload/2000123/confirm')
        response.mustcontain("forbidden")
    
    def test_no_such_experiment_should_invoke_forbidden(self):
        response = self.ptmscoutapp.get('/experiments/2000123')
        response.mustcontain("forbidden")
        
    def test_no_such_protein_should_give_info(self):
        response = self.ptmscoutapp.get('/proteins/20131234234')
        response.mustcontain(strings.error_resource_not_found_page_title)