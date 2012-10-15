from webtest.app import TestApp
from ptmscout import main
from steps.bot import Bot

def before_all(context, feature):
    settings = { 'sqlalchemy.url': "mysql+mysqldb://ptmscout_web:ptmscout1@localhost:3306/ptmscout_dev" }
    app = main({}, **settings)
    context.ptmscoutapp = TestApp(app)

def before_scenario(context, feature):
    context.active_user = Bot(context.ptmscoutapp)
    context.active_user.register()
    context.active_user.activate()

def after_scenario(context, feature):
    from ptmscout.database import DBSession
    DBSession.rollback()

def after_all(context, feature):
    from ptmscout.database import DBSession
    del context.ptmscoutapp
    DBSession.remove()
    