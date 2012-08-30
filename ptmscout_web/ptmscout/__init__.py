from pyramid.config import Configurator

from sqlalchemy import engine_from_config
from database.models import DBSession
from pyramid.authentication import AuthTktAuthenticationPolicy
from pyramid.authorization import ACLAuthorizationPolicy
from database.user import getUserByRequest

from pyramid.security import (
    Allow,
    Authenticated,
    )
from ptmscout.user_management import forbidden_view
from pyramid.exceptions import Forbidden

class RootFactory(object):
    __acl__ = [ (Allow, Authenticated, 'upload') ]
    def __init__(self, request):
        pass

def main(global_config, **settings):
    """ This function returns a Pyramid WSGI application.
    """
    config = Configurator(settings=settings,
                          root_factory='ptmscout.RootFactory')
    engine = engine_from_config(settings, 'sqlalchemy.')
    DBSession.configure(bind=engine)
    
    authn_policy = AuthTktAuthenticationPolicy('$up3r$3cr3t')
    authz_policy = ACLAuthorizationPolicy()
    
    config.set_request_property(getUserByRequest, 'user', reify=True)
    config.set_authentication_policy(authn_policy)
    config.set_authorization_policy(authz_policy)
    
    config.add_static_view('static', 'static', cache_max_age=3600)
    config.add_route('login', '/login')
    config.add_route('logout', '/logout')
    config.add_route('upload', '/upload')
    config.add_route('process_login', '/process_login')
    config.add_route('register', '/register')
    config.add_route('process_registration', '/process_registration')
    config.add_view(forbidden_view, context=Forbidden)
    
    config.scan()
    return config.make_wsgi_app()