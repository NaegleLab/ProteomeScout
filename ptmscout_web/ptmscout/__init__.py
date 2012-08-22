from pyramid.config import Configurator

from sqlalchemy import engine_from_config
from database.models import DBSession

def main(global_config, **settings):
    """ This function returns a Pyramid WSGI application.
    """
    config = Configurator(settings=settings)
    engine = engine_from_config(settings, 'sqlalchemy.')
    DBSession.configure(bind=engine)
    config.add_static_view('static', 'static', cache_max_age=3600)
    config.add_route('login', '/login')
    config.add_route('process_login', '/process_login')
    config.add_route('register', '/register')
    config.add_route('process_registration', '/process_registration')
    config.scan()
    return config.make_wsgi_app()
