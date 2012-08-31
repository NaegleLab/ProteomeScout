from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config
import urllib
import utils.webutils as webutils
from database import user
from utils import crypto
from database.user import NoSuchUser
from pyramid import security
import re
import config
from ptmscout.utils import mail, transactions



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
            'pageTitle': "User Registration"}
    
@view_config(route_name='process_registration', renderer='templates/information.pt')
def user_registration_success(request):
    result = __process_registration(request)
    if result == True:
        return {
            'pageTitle': "User Registration",
            'header': "Registration Successful",
            'message': "A confirmation e-mail has been sent to the specified e-mail address. Please check your e-mail to complete your registration."}
    else:
        raise HTTPFound(request.application_url+"/register?"+urllib.urlencode(result))



@view_config(route_name='login', renderer='templates/login.pt')
def user_login(request):
    username = webutils.get(request, 'username', "")
    reason = webutils.get(request, 'reason', None)
    
    return {
            'username': username,
            'reason':reason,
            'pageTitle': "Login"}


@view_config(route_name='process_login', renderer='templates/information.pt')
def user_login_success(request):
    result = __process_login(request)
    if result == True:
        return {
                'pageTitle': "Login",
                'header': "Login Successful",
                'message': "You have successfully logged in."}
    else:
        raise HTTPFound(request.application_url+"/login?"+urllib.urlencode(result))

@view_config(route_name='logout', renderer='templates/information.pt')
def user_logout(request):    
    request.response.headers.extend(security.forget(request))
    request.user = None
    
    return {
            'pageTitle': "Logout",
            'header': "Logout Successful",
            'message': "You have successfully logged out."}
    
    
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
    
    transactions.commit()    
    
    return {'pageTitle':"Account Activation",
            'header':"Account Activation Succeeded",
            'message':"Your account is now active. Please <a href=\""+request.application_url+"/login\">login</a>"}

## Internal Functions

def __failed_account_activation(request):
    return {'pageTitle':"Account Activation",
            'header':"Account Activation Failure",
            'message':"The specified account is not valid, please try <a href=\""+request.application_url+"/register\">registering</a>"}

def __process_login(request):
    username = webutils.post(request, 'username', "")
    password = webutils.post(request, 'password', "")
    
    resp_dict = {'username': username}
    
    if username == "" or password == "":
        resp_dict['reason'] = "All fields are required"
        return resp_dict
    
    try:
        ptm_user = user.getUserByUsername(username)
        
        _, salted_password = crypto.saltedPassword(password, ptm_user.salt)
        
        if salted_password != ptm_user.salted_password:
            raise NoSuchUser()
        
        if not ptm_user.active:
            resp_dict['reason'] = "Account has not been activated"
            return resp_dict
        
        request.user = ptm_user
        request.response.headers.extend(security.remember(request, ptm_user.username))
    except NoSuchUser:
        resp_dict['reason'] = "Credentials incorrect"
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
        resp_dict['reason'] = "Form fields cannot be empty"
        return resp_dict

    try:
        _ptm_user = user.getUserByUsername(username)
        resp_dict['reason'] = "Username is already in use"
        return resp_dict
    except NoSuchUser:
        pass
    
    
    email_regex = re.compile("[a-z0-9\.\-\_]+@[a-z0-9\.\-\_]+\.([a-z]+)$", re.I)
    
    matches = email_regex.match(email)
    if not matches:
        resp_dict['reason'] = "Email address is invalid"
        return resp_dict

    domain = matches.group(1)
    if domain != "edu":
        resp_dict['reason'] = "Email address must belong to .edu domain"
        return resp_dict

    if len(pass1) < config.MINIMUM_PASSWORD_LENGTH:
        resp_dict['reason'] = "Password must be at least %d characters in length" % config.MINIMUM_PASSWORD_LENGTH
        return resp_dict

    if( pass1 != pass2 ):
        resp_dict['reason'] = "Password confirmation does not match"
        return resp_dict

    ptm_user = user.PTMUser(username, name, email, institute)
    ptm_user.createUser(pass1)
    ptm_user.saveUser()
    
    mail.send_automail_message(request, [email], "PTMScout Account Activiation Details", "%s, \n\nThank you for choosing PTMScout for your research.\n\nYou can activate your new account by visiting <a href=\"%s/activate_account?username=%s&token=%s\">this link</a>.\n\nThanks,\n-The PTMScout Team" % (ptm_user.name, request.application_url, ptm_user.username, ptm_user.activation_token))
    
    transactions.commit()

    return True