from pyramid.view import forbidden_view_config, view_config
from smtplib import SMTPRecipientsRefused
from sqlalchemy.exc import IntegrityError
from ptmscout.database.experiment import ExperimentAccessForbidden,\
    NoSuchExperiment, ExperimentNotAvailable
from ptmscout.database.protein import NoSuchProtein
from ptmscout.config import strings
from ptmscout.database.upload import NoSuchSession, SessionAccessForbidden
import logging

@view_config(context=IntegrityError, renderer='ptmscout:templates/info/information.pt')
def internal_database_error_view(exc, request):
    log = logging.getLogger('ptmscout')
    log.exception(exc)
    
    return {'pageTitle':"Database Error",
            'header':"Database Error",
            'message':"There was an internal database error during the processing of your request, please try again later.",
            'redirect': request.application_url + "/experiments"}

@view_config(context=SMTPRecipientsRefused, renderer='ptmscout:templates/info/information.pt')
def internal_mailer_error_view(exc, request):
    log = logging.getLogger('ptmscout')
    log.exception(exc)
    
    return {'pageTitle':"Mailer Error",
            'header':"Mailer Error",
            'message':"An error was encountered while processing outgoing SMTP messages, please try again later.",
            'redirect': request.application_url + "/experiments"}

@view_config(context=NoSuchProtein, renderer='ptmscout:templates/info/information.pt')
def redirect_no_such_protein(exc, request):
    return {'pageTitle':strings.error_resource_not_found_page_title,
            'header':strings.error_resource_not_found_page_title,
            'message':strings.error_protein_not_found_message,
            'redirect': request.application_url + "/proteins"}

@view_config(context=NoSuchSession, renderer='ptmscout:templates/info/forbidden.pt')
def redirect_no_such_session(exc, request):
    return {'pageTitle': "Forbidden"}

@view_config(context=SessionAccessForbidden, renderer='ptmscout:templates/info/forbidden.pt')
def redirect_session_forbidden(exc, request):
    return {'pageTitle': "Forbidden"}

@view_config(context=NoSuchExperiment, renderer='ptmscout:templates/info/forbidden.pt')
def redirect_no_such_experiment(exc, request):
    return {'pageTitle': "Forbidden"}

@view_config(context=ExperimentNotAvailable, renderer='ptmscout:templates/info/information.pt')
def experiment_not_available(exc, request):
    return {'pageTitle':strings.error_resource_not_ready_page_title,
            'header':strings.error_resource_not_ready_page_title,
            'message':strings.error_resource_not_ready_message % (request.application_url + "/account/experiments")}

@view_config(context=ExperimentAccessForbidden, renderer='ptmscout:templates/info/forbidden.pt')
def redirect_forbidden(request):
    return {'pageTitle': "Forbidden"}

@forbidden_view_config(renderer='ptmscout:templates/info/forbidden.pt')
def forbidden_view(request):
    return {'pageTitle': "Forbidden"}