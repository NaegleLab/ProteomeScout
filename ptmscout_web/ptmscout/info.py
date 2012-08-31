from pyramid.view import view_config

@view_config(route_name='about', renderer='templates/about.pt')
def about_view(request):
    return {'pageTitle':"About PTMScout"}

@view_config(route_name='terms', renderer='templates/terms.pt')
def terms_view(request):
    return {'pageTitle':"PTMScout Terms of Use"}