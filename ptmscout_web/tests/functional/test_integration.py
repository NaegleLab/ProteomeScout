import unittest
from webtest.app import TestApp
from ptmscout import main

class InfoFunctionalTests(unittest.TestCase):
    def setUp(self):
        settings = { 'sqlalchemy.url': "mysql+mysqldb://ptmscout_web:ptmscout1@localhost:3306/ptmscout_dev" }
        app = main({}, **settings)
        self.ptmscoutapp = TestApp(app)
        
    def tearDown(self):
        from ptmscout.database import DBSession
        del self.ptmscoutapp
        DBSession.remove()

    def test_about_renderer(self):
        self.ptmscoutapp.get('/about', status=200)
    
    def test_terms_renderer(self):
        self.ptmscoutapp.get('/terms', status=200)
        
    def test_integrated_protein_views(self):
        self.ptmscoutapp.get('/proteins/35546/modifications')
        self.ptmscoutapp.get('/proteins/35546/expression')
        self.ptmscoutapp.get('/proteins/35546/GO')
        
    