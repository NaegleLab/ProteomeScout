from pyramid.view import forbidden_view_config, view_config
from ptmscout.database import InternalDatabaseError


@view_config(context=InternalDatabaseError, renderer='templates/information.pt')
def internal_database_error_view(exc, request):
    return {'pageTitle':"Database Error",
            'header':"Database Error",
            'message':"There was an internal error during the processing of your request, please try again later."}

@forbidden_view_config(renderer='templates/forbidden.pt')
def forbidden_view(request):
    return {'pageTitle': "Forbidden"}