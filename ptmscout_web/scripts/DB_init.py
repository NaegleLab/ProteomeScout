from sqlalchemy import engine_from_config
from paste.deploy.loadwsgi import appconfig
from ptmscout.database import DBSession, Base  # base declarative object
from pyramid import testing
import sys, logging

class DatabaseInitialization():
    @classmethod
    def setUpClass(cls, settings_path, sql_echo=False):
        settings = appconfig('config:%s' % (settings_path))
        cls.engine = engine_from_config(settings, prefix='sqlalchemy.')
        settings['sqlalchemy.echo'] = sql_echo
        settings['sqlalchemy.echo_pool'] = sql_echo

    def log_capture(self, log_stream, level=logging.DEBUG):
        logger = logging.getLogger(log_stream)
        logger.level = level
        stream_handler = logging.StreamHandler(sys.stdout)
        logger.addHandler(stream_handler)

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
