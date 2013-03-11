import unittest
from pyramid import testing
from paste.deploy.loadwsgi import appconfig

from sqlalchemy import engine_from_config
from ptmscout.database import DBSession, Base  # base declarative object

import os
settings = appconfig('config:/' + os.path.join('data','ptmscout','ptmscout_web', 'test.ini'))

class DBTestCase(unittest.TestCase):
    engine=None

    @classmethod
    def setUpClass(cls):
        if cls.engine == None:
            cls.engine = engine_from_config(settings, prefix='sqlalchemy.')

    def setUp(self):
        connection = self.engine.connect()

        # begin a non-ORM transaction
        self.trans = connection.begin()

        # bind an individual Session to the connection
        DBSession.configure(bind=connection)
        self.session = DBSession
        Base.session = self.session

    def tearDown(self):
        # rollback - everything that happened with the
        # Session above (including calls to commit())
        # is rolled back.
        testing.tearDown()
        self.trans.rollback()
        self.session.close()
