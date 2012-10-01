from webtest.app import TestApp
from ptmscout import main
from steps.bot import Bot

def before_feature(context, feature):
    settings = { 'sqlalchemy.url': "mysql+mysqldb://ptmscout_web:ptmscout1@localhost:3306/ptmscout_dev" }
    print "Starting test app"
    app = main({}, **settings)
    context.ptmscoutapp = TestApp(app)

def before_scenario(context, feature):
    context.active_user = Bot(context.ptmscoutapp)
    context.active_user.register()
    context.active_user.activate()

def after_scenario(context, feature):
    from ptmscout.database import DBSession
    DBSession.rollback()

def after_feature(context, feature):
    from ptmscout.database import DBSession
    del context.ptmscoutapp
    DBSession.remove()
    