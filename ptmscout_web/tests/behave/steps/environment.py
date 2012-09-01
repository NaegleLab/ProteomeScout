from webtest.app import TestApp
from ptmscout import main

def before_feature(context, feature):
    settings = { 'sqlalchemy.url': "mysql+mysqldb://ptmscout_web:ptmscout1@localhost:3306/ptmscout" }
    print "Starting test app"
    app = main({}, **settings)
    context.ptmscoutapp = TestApp(app)


def after_feature(context, feature):
    from ptmscout.database import DBSession
    del context.ptmscoutapp
    DBSession.remove()