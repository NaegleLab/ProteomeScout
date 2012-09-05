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
from ptmscout.database import DBSession

@view_config(route_name='account_management', renderer='templates/account.pt', permission='private')
def manage_account(request):
    reason = webutils.get(request, 'reason', None)
    return {'username': request.user.username,
            'fullname': request.user.name,
            'email': request.user.email,
            'institution': request.user.institution,
            'pageTitle': "Account Management",
            'reason': reason}

@view_config(route_name='change_password', permission='private')
def change_password(request):
    oldpass = webutils.post(request, 'old_pass', "")
    newpass1 = webutils.post(request, 'new_pass1', "")
    newpass2 = webutils.post(request, 'new_pass2', "")
    
    if oldpass == "" or newpass1 == "" or newpass2 == "":
        raise HTTPFound(request.application_url + "/account?" + urllib.urlencode({'reason':"Form fields cannot be empty"}))
    
    if newpass1 != newpass2:
        raise HTTPFound(request.application_url + "/account?" + urllib.urlencode({'reason':"New password confirmation did not match"}))
    
    _, old_salted_pass = crypto.saltedPassword(oldpass, request.user.salt)
    
    if old_salted_pass != request.user.salted_password:
        raise HTTPFound(request.application_url + "/account?" + urllib.urlencode({'reason':"Supplied password was incorrect"}))
    
    _, new_salted_pass = crypto.saltedPassword(newpass1, request.user.salt)
    request.user.salted_password = new_salted_pass
    request.user.saveUser()
    
    transactions.commit()
    
    raise HTTPFound(request.application_url + "/change_password_success")
    
@view_config(route_name='change_password_success', renderer='templates/information.pt', permission='private')
def change_password_success(request):
    return {'pageTitle': "Change Password",
            'message': "Password successfully changed.",
            'header': "Success"}

@view_config(route_name='forgot_password', renderer='templates/forgot_password.pt')
def forgot_password(request):
    email = webutils.get(request, 'email', "")
    reason = webutils.get(request, 'reason', None)
    
    return {'pageTitle': "Forgotten Password Retrieval",
            'email': email,
            'reason': reason}

@view_config(route_name='process_forgot_password', renderer='templates/information.pt')
def process_forgot_password(request):
    email = webutils.post(request, 'email', "")
    
    if email == "":
        raise HTTPFound(request.application_url + "/forgot_password?" + urllib.urlencode({'email':email, 'reason':"Form fields cannot be empty"}))
    
    new_password = crypto.randomString(5)
    
    try:
        ptm_user = user.getUserByEmail(email)
        _, ptm_user.salted_password = crypto.saltedPassword(new_password, ptm_user.salt)
        ptm_user.saveUser()
    except NoSuchUser:
        raise HTTPFound(request.application_url + "/forgot_password?" + urllib.urlencode({'email':email, 'reason': "E-mail address does not match any user record"}))
    
    login_url = request.application_url + "/login"
    account_url = request.application_url + "/account"
    
    message = """%s,
        
        Your password in PTMScout has been reset, your new login credentials are:
        Username: %s
        Password: %s
        
        Please visit <a href="%s">PTMScout</a> to login.
        After logging in, your can change your password <a href="%s">here</a>.
        
        -PTMScout Administrator
        """ % (ptm_user.name, ptm_user.username, new_password, login_url, account_url)
    
    mail.send_automail_message(request, [email], "PTMScout password reset", message)
    
    transactions.commit()
    
    return {'pageTitle': "Forgotten Password Retrieval",
            'header': "Password Reset Success",
            'message':"Your username and a temporary password have been sent to your e-mail address"}

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

    ptm_user = user.User(username, name, email, institute)
    ptm_user.createUser(pass1)
    ptm_user.saveUser()
    
    mail.send_automail_message(request, [email], "PTMScout Account Activiation Details", "%s, \n\nThank you for choosing PTMScout for your research.\n\nYou can activate your new account by visiting <a href=\"%s/activate_account?username=%s&token=%s\">this link</a>.\n\nThanks,\n-The PTMScout Team" % (ptm_user.name, request.application_url, ptm_user.username, ptm_user.activation_token))
    
    transactions.commit()

    return True