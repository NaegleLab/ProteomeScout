import unittest
from pyramid import testing
from ptmscout import main
from webtest.app import TestApp
from tests.behave.steps.bot import Bot
from paste.deploy.loadwsgi import appconfig
import os, sys
import logging

class IntegrationTestCase(unittest.TestCase):
    def setUp(self):
        settings = appconfig('config:/' + os.path.join('data','ptmscout','ptmscout_web', 'test.ini'))
        app = main({}, **settings)
        self.ptmscoutapp = TestApp(app)
        self.bot = Bot(self.ptmscoutapp)
        self.bot.register_and_login()

    def log_capture(self):
        logger = logging.getLogger('ptmscout')
        logger.level = logging.DEBUG
        stream_handler = logging.StreamHandler(sys.stdout)
        logger.addHandler(stream_handler)
       
    def tearDown(self):
        from ptmscout.database import DBSession
        del self.ptmscoutapp
        DBSession.remove()
        
class UnitTestCase(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        
    def tearDown(self):
        testing.tearDown()
        
