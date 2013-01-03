from sqlalchemy import engine_from_config
from ptmscout.database import DBSession, Base  # base declarative object
from pyramid import testing

class DatabaseInitialization():
    @classmethod
    def setUpClass(cls, settings):
        cls.engine = engine_from_config(settings, prefix='sqlalchemy.')

    def setUp(self):
        self.connection = self.engine.connect()

        # begin a non-ORM transaction
        self.trans = self.connection.begin()

        # bind an individual Session to the connection
        DBSession.configure(bind=self.connection)
        self.session = DBSession
        Base.session = self.session

    def new_transaction(self):
        self.trans = self.connection.begin()

    def commit(self):
        self.trans.commit()

    def rollback(self):
        testing.tearDown()
        self.trans.rollback()
        self.session.close()

    def tearDown(self):
        testing.tearDown()
        self.trans.commit()
        self.session.close()
