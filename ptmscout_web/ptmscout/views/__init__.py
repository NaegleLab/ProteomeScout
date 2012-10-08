from pyramid.exceptions import Forbidden
from ptmscout.views.errors import forbidden_view

def add_views(config):
    config.include('pyramid_mailer')
    
    config.add_static_view('static', 'static', cache_max_age=3600)
    
    config.add_route('about', '/about')
    config.add_route('terms', '/terms')
    
    config.add_route('redirect_to_experiments','/')
    config.add_route('experiments','/experiments')
    config.add_route('experiment','/experiments/{id}')
    config.add_route('experiment_GO','/experiments/{id}/GO')
    config.add_route('experiment_pfam','/experiments/{id}/pfam')
    config.add_route('experiment_predictions','/experiments/{id}/predictions')
    config.add_route('experiment_summary','/experiments/{id}/summary')
    config.add_route('experiment_browse','/experiments/{id}/browse')
    config.add_route('upload', '/upload')
    
    config.add_route('protein_data', '/proteins/{id}/data')
    config.add_route('protein_GO', '/proteins/{id}/GO')
    config.add_route('protein_expression', '/proteins/{id}/expression')
    config.add_route('protein_mod_sites', '/proteins/{id}/modifications')
    config.add_route('protein_main', '/proteins/{id}')
    config.add_route('protein_search', '/proteins')
    
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
    config.add_route('invite_experiment', '/account/experiments/{id}/invite')
    config.add_route('share_experiment', '/account/experiments/{id}/share')
    config.add_route('publish_experiment', '/account/experiments/{id}/publish')
    config.add_route('privatize_experiment', '/account/experiments/{id}/unpublish')
    
    config.add_route('change_password', '/change_password')
    config.add_route('change_password_success', '/change_password_success')
    
    config.add_view(forbidden_view, context=Forbidden)