import unittest
from pyramid import testing
from ptmscout import main
from webtest.app import TestApp
from tests.behave.steps.bot import Bot
from paste.deploy.loadwsgi import appconfig
import os

class IntegrationTestCase(unittest.TestCase):
    def setUp(self):
        settings = appconfig('config:/' + os.path.join('data','ptmscout','ptmscout_web', 'test.ini'))
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
        