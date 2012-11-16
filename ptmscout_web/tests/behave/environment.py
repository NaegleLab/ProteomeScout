from webtest.app import TestApp
from ptmscout import main
from steps.bot import Bot
from paste.deploy.loadwsgi import appconfig
import os

def before_all(context):
    settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'test.ini'))
    
    app = main({}, **settings)
    context.ptmscoutapp = TestApp(app)
    
def before_scenario(context, feature):
    context.active_user = Bot(context.ptmscoutapp)
    context.active_user.register()
    context.active_user.activate()

def after_scenario(context, feature):
    from ptmscout.database import DBSession
    DBSession.rollback()

def after_all(context):
    from ptmscout.database import DBSession
    del context.ptmscoutapp
    DBSession.remove()
    