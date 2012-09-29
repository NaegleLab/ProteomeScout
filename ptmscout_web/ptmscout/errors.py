from pyramid.view import forbidden_view_config, view_config
from smtplib import SMTPRecipientsRefused
from sqlalchemy.exc import IntegrityError
from ptmscout.database.experiment import ExperimentAccessForbidden,\
    NoSuchExperiment
from pyramid.httpexceptions import HTTPForbidden

@view_config(context=IntegrityError, renderer='templates/information.pt')
def internal_database_error_view(exc, request):
    return {'pageTitle':"Database Error",
            'header':"Database Error",
            'message':"There was an internal database error during the processing of your request, please try again later.",
            'redirect': request.application_url + "/experiments"}

@view_config(context=SMTPRecipientsRefused, renderer='templates/information.pt')
def internal_mailer_error_view(exc, request):
    return {'pageTitle':"Mailer Error",
            'header':"Mailer Error",
            'message':"An error was encountered while processing outgoing SMTP messages, please try again later.",
            'redirect': request.application_url + "/experiments"}

@view_config(context=NoSuchExperiment, renderer='templates/forbidden.pt')
def redirect_no_such_experiment(exc, request):
    return {'pageTitle': "Forbidden"}

@view_config(context=ExperimentAccessForbidden, renderer='templates/forbidden.pt')
def redirect_forbidden(request):
    return {'pageTitle': "Forbidden"}

@forbidden_view_config(renderer='templates/forbidden.pt')
def forbidden_view(request):
    return {'pageTitle': "Forbidden"}