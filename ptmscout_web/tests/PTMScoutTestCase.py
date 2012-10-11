import unittest
from pyramid import testing
from ptmscout import main
from webtest.app import TestApp
from tests.behave.steps.bot import Bot

class IntegrationTestCase(unittest.TestCase):
    def setUp(self):
        settings = { 'sqlalchemy.url': "mysql+mysqldb://ptmscout_web:ptmscout1@localhost:3306/ptmscout_dev" }
        app = main({}, **settings)
        self.ptmscoutapp = TestApp(app)
        self.bot = Bot(self.ptmscoutapp)
        self.bot.register_and_login()
        
    def tearDown(self):
        from ptmscout.database import DBSession
        del self.ptmscoutapp
        DBSession.remove()
        
class UnitTestCase(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        
    def tearDown(self):
        testing.tearDown()
        