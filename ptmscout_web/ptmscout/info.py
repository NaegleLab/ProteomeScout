from pyramid.view import view_config
from ptmscout import strings

@view_config(route_name='about', renderer='templates/about.pt')
def about_view(request):
    return {'pageTitle':strings.about_page_title}

@view_config(route_name='terms', renderer='templates/terms.pt')
def terms_view(request):
    return {'pageTitle':strings.view_page_title}