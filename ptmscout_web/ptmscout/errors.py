from pyramid.view import forbidden_view_config, view_config
from ptmscout.utils.transactions import InternalDatabaseError, MailerError


@view_config(context=InternalDatabaseError, renderer='templates/information.pt')
def internal_database_error_view(exc, request):
    return {'pageTitle':"Database Error",
            'header':"Database Error",
            'message':"There was an internal error during the processing of your request, please try again later."}

@view_config(context=MailerError, renderer='templates/information.pt')
def internal_mailer_error_view(exc, request):
    return {'pageTitle':"Mailer Error",
            'header':"Mailer Error",
            'message':"An error was encountered while processing outgoing SMTP messages, please try again later."}


@forbidden_view_config(renderer='templates/forbidden.pt')
def forbidden_view(request):
    return {'pageTitle': "Forbidden"}