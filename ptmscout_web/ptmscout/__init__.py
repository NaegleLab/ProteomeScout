from pyramid.config import Configurator
from sqlalchemy import engine_from_config
from database import DBSession
from pyramid.authentication import AuthTktAuthenticationPolicy
from pyramid.authorization import ACLAuthorizationPolicy
from database.user import getUserByRequest
import views

from pyramid.security import (
    Allow,
    Authenticated,
    )
from ptmscout.views.errors import forbidden_view
from pyramid.exceptions import Forbidden
from ptmscout.database.user import getUsernameByRequest
from ptmscout import webservice

engine = None

class RootFactory(object):
    __acl__ = [ (Allow, Authenticated, 'private') ]
    def __init__(self, request):
        pass
    
def includeme(config):
    views.add_views(config)

def main(global_config, **settings):
    global engine
    """ This function returns a Pyramid WSGI application.
    """
    config = Configurator(settings=settings,
                          root_factory='ptmscout.RootFactory')

    if engine == None:
        engine = engine_from_config(settings, 'sqlalchemy.')
        DBSession.configure(bind=engine)
    
    authn_policy = AuthTktAuthenticationPolicy('$up3r$3cr3t')
    authz_policy = ACLAuthorizationPolicy()
    
    config.set_request_property(getUserByRequest, 'user', reify=True)
    config.set_request_property(getUsernameByRequest, 'username', reify=True)
    config.set_authentication_policy(authn_policy)
    config.set_authorization_policy(authz_policy)
    
    views.add_views(config)
    webservice.add_views(config)
    
    config.scan()
    return config.make_wsgi_app()
