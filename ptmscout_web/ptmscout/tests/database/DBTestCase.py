import unittest
from pyramid import testing
from paste.deploy.loadwsgi import appconfig

from sqlalchemy import engine_from_config
from sqlalchemy.orm import sessionmaker
from database.models import DBSession, Base  # base declarative object

import os
here = os.path.dirname(__file__)
settings = appconfig('config:' + os.path.join(here, '../../../', 'development.ini'))

class DBTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = engine_from_config(settings, prefix='sqlalchemy.')
        cls.Session = sessionmaker()

    def setUp(self):
        connection = self.engine.connect()

        # begin a non-ORM transaction
        self.trans = connection.begin()

        # bind an individual Session to the connection
        DBSession.configure(bind=connection)
        self.session = self.Session(bind=connection)
        Base.session = self.session

    def tearDown(self):
        # rollback - everything that happened with the
        # Session above (including calls to commit())
        # is rolled back.
        testing.tearDown()
        self.trans.rollback()
        self.session.close()