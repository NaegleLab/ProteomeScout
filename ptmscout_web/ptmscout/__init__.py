from pyramid.config import Configurator

from sqlalchemy import engine_from_config
from database import DBSession
from pyramid.authentication import AuthTktAuthenticationPolicy
from pyramid.authorization import ACLAuthorizationPolicy
from database.user import getUserByRequest

from pyramid.security import (
    Allow,
    Authenticated,
    )
from ptmscout.errors import forbidden_view
from pyramid.exceptions import Forbidden

class RootFactory(object):
    __acl__ = [ (Allow, Authenticated, 'private') ]
    def __init__(self, request):
        pass
    
def add_views(config):
    config.include('pyramid_mailer')
    
    config.add_static_view('static', 'static', cache_max_age=3600)
    
    config.add_route('about', '/about')
    config.add_route('terms', '/terms')
    
    config.add_route('redirect_to_experiments','/')
    config.add_route('experiments','/experiments')
    config.add_route('experiment','/experiments/{id}')
    config.add_route('upload', '/upload')
    
    config.add_route('login', '/login')
    config.add_route('process_login', '/process_login')
    config.add_route('logout', '/logout')
    
    config.add_route('register', '/register')
    config.add_route('process_registration', '/process_registration')
    config.add_route('activate_account', '/activate_account')
    
    config.add_route('forgot_password', '/forgot_password')
    config.add_route('process_forgot_password', '/retrieve_password')
    
    config.add_route('account_management', '/account')
    config.add_route('my_experiments', '/account/experiments')
    config.add_route('share_experiment', '/account/experiments/{id}/share')
#    config.add_route('publish_experiment', '/account/experiments/{id}/publish')
    
    config.add_route('change_password', '/change_password')
    config.add_route('change_password_success', '/change_password_success')
    
    config.add_view(forbidden_view, context=Forbidden)
    
def includeme(config):
    add_views(config)

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
    
    add_views(config)
    
    config.scan()
    return config.make_wsgi_app()