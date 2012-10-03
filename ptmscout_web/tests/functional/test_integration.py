import unittest
from webtest.app import TestApp
from ptmscout import main
from ptmscout.config import strings
from tests.PTMScoutTestCase import IntegrationTestCase

class InfoFunctionalTests(IntegrationTestCase):
    def test_about_renderer(self):
        self.ptmscoutapp.get('/about', status=200)
    
    def test_terms_renderer(self):
        self.ptmscoutapp.get('/terms', status=200)
        
    def test_integrated_protein_views(self):
        self.ptmscoutapp.get('/proteins')
        self.ptmscoutapp.get('/proteins/35546/modifications')
        self.ptmscoutapp.get('/proteins/35546/expression')
        self.ptmscoutapp.get('/proteins/35546/GO')
        
    def test_integrated_experiment_views(self):
        self.ptmscoutapp.get('/experiments', status=200)
        self.ptmscoutapp.get('/experiments/26')
        self.ptmscoutapp.get('/experiments/26/summary')
        self.ptmscoutapp.get('/experiments/26/browse')
        self.ptmscoutapp.post('/experiments/26/browse', {'submitted':"true", 'acc_search':"ACK1", 'stringency':"1"})
        
    def test_login_views(self):
        self.ptmscoutapp.get('/login')
        self.ptmscoutapp.get('/register')
        
    def test_forbidden_should_invoke_on_unauthorized_access(self):
        response = self.ptmscoutapp.get('/upload')
        response.mustcontain("forbidden")
        
    def test_no_such_experiment_should_invoke_forbidden(self):
        response = self.ptmscoutapp.get('/experiments/2000123')
        response.mustcontain("forbidden")
        
    def test_no_such_protein_should_give_info(self):
        response = self.ptmscoutapp.get('/proteins/20131234234')
        
        response.mustcontain(strings.error_resource_not_found_page_title)