from pyramid.view import view_config
from ptmscout.config import strings

@view_config(route_name='about', renderer='ptmscout:templates/info/about.pt')
def about_view(request):
    return {'pageTitle':strings.about_page_title}

@view_config(route_name='terms', renderer='ptmscout:templates/info/terms.pt')
def terms_view(request):
    return {'pageTitle':strings.view_page_title}

@view_config(route_name='portal', renderer='ptmscout:templates/front_page.pt')
def portal_view(request):
    return {'pageTitle':strings.portal_page_title}
