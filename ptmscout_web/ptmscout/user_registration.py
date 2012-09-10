from database import user
from database.user import NoSuchUser
from ptmscout.utils import mail
from pyramid import security
from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config
from utils import crypto
import config
import re
import strings
import urllib
import utils.webutils as webutils

@view_config(route_name='forgot_password', renderer='templates/forgot_password.pt')
def forgot_password(request):
    email = webutils.get(request, 'email', "")
    reason = webutils.get(request, 'reason', None)
    
    return {'pageTitle': strings.forgotten_password_page_title,
            'email': email,
            'reason': reason}

@view_config(route_name='process_forgot_password', renderer='templates/information.pt')
def process_forgot_password(request):
    email = webutils.post(request, 'email', "")
    
    if email == "":
        raise HTTPFound(request.application_url + "/forgot_password?" + urllib.urlencode({'email':email, 'reason':strings.failure_reason_form_fields_cannot_be_empty}))
    
    new_password = crypto.randomString(5)
    
    try:
        ptm_user = user.getUserByEmail(email)
        _, ptm_user.salted_password = crypto.saltedPassword(new_password, ptm_user.salt)
        ptm_user.saveUser()
    except NoSuchUser:
        raise HTTPFound(request.application_url + "/forgot_password?" + urllib.urlencode({'email':email, 'reason': strings.failure_reason_email_address_not_on_record}))
    
    login_url = request.application_url + "/login"
    account_url = request.application_url + "/account"
    
    message = strings.forgotten_password_email_message % (ptm_user.name, ptm_user.username, new_password, login_url, account_url)
    
    mail.send_automail_message(request, [email], strings.forgotten_password_email_subject, message)
    
    return {'pageTitle': strings.forgotten_password_page_title,
            'header': strings.forgotten_password_success_header,
            'message': strings.forgotten_password_success_message}


@view_config(route_name='register',renderer='templates/user_registration.pt')
def user_registration_view(request):
    username = webutils.get(request, 'username', "")  
    email = webutils.get(request, 'email', "")
    name = webutils.get(request, 'name', "")
    institution = webutils.get(request, 'institution', "")
    reason = webutils.get(request, 'reason', None)
    
    return {
            'username':username,
            'email':email,
            'name':name,
            'institution':institution,
            'reason':reason,
            'pageTitle': strings.user_registration_page_title}
    
@view_config(route_name='process_registration', renderer='templates/information.pt')
def user_registration_success(request):
    result = __process_registration(request)
    if result == True:
        return {
            'pageTitle': strings.user_registration_page_title,
            'header': strings.user_registration_success_header,
            'message': strings.user_registration_success_message}
    else:
        raise HTTPFound(request.application_url+"/register?"+urllib.urlencode(result))



@view_config(route_name='login', renderer='templates/login.pt')
def user_login(request):
    username = webutils.get(request, 'username', "")
    reason = webutils.get(request, 'reason', None)
    
    return {
            'username': username,
            'reason':reason,
            'pageTitle': strings.login_page_title}


@view_config(route_name='process_login', renderer='templates/information.pt')
def user_login_success(request):
    result = __process_login(request)
    if result == True:
        return {
                'pageTitle': strings.login_page_title,
                'header': strings.login_page_success_header,
                'message': strings.login_page_success_message}
    else:
        raise HTTPFound(request.application_url+"/login?"+urllib.urlencode(result))

@view_config(route_name='logout', renderer='templates/information.pt')
def user_logout(request):    
    request.response.headers.extend(security.forget(request))
    request.user = None
    
    return {
            'pageTitle': strings.logout_page_title,
            'header': strings.logout_page_header,
            'message': strings.logout_page_message}
    
    
@view_config(route_name='activate_account', renderer='templates/information.pt')
def user_account_activation(request):
    username = webutils.get(request, 'username', "").strip()
    token = webutils.get(request, 'token', "").strip()
    
    if username == "" or token == "":
        return __failed_account_activation(request)
    
    ptm_user = None
    try:
        ptm_user = user.getUserByUsername(username)
    except NoSuchUser:
        return __failed_account_activation(request)
    
    if ptm_user.activation_token != token:
        return __failed_account_activation(request)
    
    ptm_user.setActive()
    ptm_user.saveUser()
    
    return {'pageTitle': strings.account_activation_page_title,
            'header': strings.account_activation_success_header,
            'message': strings.account_activation_success_message % (request.application_url+"/login")}

## Internal Functions

def __failed_account_activation(request):
    return {'pageTitle':strings.account_activation_page_title,
            'header':strings.account_activation_failed_header,
            'message':strings.account_activation_failed_message % (request.application_url+"/register")}

def __process_login(request):
    username = webutils.post(request, 'username', "")
    password = webutils.post(request, 'password', "")
    
    resp_dict = {'username': username}
    
    if username == "" or password == "":
        resp_dict['reason'] = strings.failure_reason_form_fields_cannot_be_empty
        return resp_dict
    
    try:
        ptm_user = user.getUserByUsername(username)
        
        _, salted_password = crypto.saltedPassword(password, ptm_user.salt)
        
        if salted_password != ptm_user.salted_password:
            raise NoSuchUser()
        
        if not ptm_user.active:
            resp_dict['reason'] = strings.failure_reason_inactive_account
            return resp_dict
        
        request.user = ptm_user
        request.response.headers.extend(security.remember(request, ptm_user.username))
    except NoSuchUser:
        resp_dict['reason'] = strings.failure_reason_incorrect_credentials
        return resp_dict
    
    return True

def __process_registration(request):
    username = webutils.post(request, 'username', "").strip()
    pass1    = webutils.post(request, 'pass1', "").strip()
    pass2    = webutils.post(request, 'pass2', "").strip()
    
    name      = webutils.post(request, 'name', "").strip()
    email     = webutils.post(request, 'email', "").strip()
    institute = webutils.post(request, 'institution', "").strip()
    
    resp_dict = {'username': username, 'email':email, 'name':name, 'institution':institute}
    
    if username == "" or pass1 == "" or pass2 == "" or name == "" or email == "" or institute == "":
        resp_dict['reason'] = strings.failure_reason_form_fields_cannot_be_empty
        return resp_dict

    try:
        _ptm_user = user.getUserByUsername(username)
        resp_dict['reason'] = strings.failure_reason_username_inuse
        return resp_dict
    except NoSuchUser:
        pass
    
    
    email_regex = re.compile(strings.email_regex, re.I)
    
    matches = email_regex.match(email)
    if not matches:
        resp_dict['reason'] = strings.failure_reason_email_not_valid
        return resp_dict

    domain = matches.group(1)
    if domain != "edu":
        resp_dict['reason'] = strings.failure_reason_email_not_academic
        return resp_dict

    if len(pass1) < config.MINIMUM_PASSWORD_LENGTH:
        resp_dict['reason'] = strings.failure_reason_password_too_short % config.MINIMUM_PASSWORD_LENGTH
        return resp_dict

    if( pass1 != pass2 ):
        resp_dict['reason'] = strings.failure_reason_new_passwords_not_matching
        return resp_dict

    ptm_user = user.User(username, name, email, institute)
    ptm_user.createUser(pass1)
    ptm_user.saveUser()
    
    mail.send_automail_message(request, 
    [email], 
    strings.user_registration_email_subject, 
    strings.user_registration_email_message % (ptm_user.name, request.application_url, ptm_user.username, ptm_user.activation_token))
    
    return True