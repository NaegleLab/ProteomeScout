from layout import site_layout
from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config, forbidden_view_config
import urllib
import utils.webutils as webutils
from database import user
from utils import crypto
from database.user import NoSuchUser
from pyramid import security



@view_config(route_name='register',renderer='templates/user_registration.pt')
def user_registration_view(request):
    username = webutils.get(request, 'username', "")  
    email = webutils.get(request, 'email', "")
    name = webutils.get(request, 'name', "")
    institution = webutils.get(request, 'institution', "")
    reason = webutils.get(request, 'reason', None)
    
    return {'layout': site_layout(),
            'username':username,
            'email':email,
            'name':name,
            'institution':institution,
            'reason':reason,
            'pageTitle': "User Registration"}
    
@forbidden_view_config(renderer='templates/forbidden.pt')
def forbidden_view(request):
    return {'layout': site_layout(),
            'pageTitle': "Forbidden"}

@view_config(route_name='process_registration', renderer='templates/information.pt')
def user_registration_success(request):
    result = __process_registration(request)
    if result == True:
        return {'layout': site_layout(),
            'pageTitle': "User Registration",
            'header': "Registration Successful",
            'message': "A confirmation e-mail has been sent to the specified e-mail address. Please check your e-mail to complete your registration."}
    else:
        raise HTTPFound(request.application_url+"/register?"+urllib.urlencode(result))



@view_config(route_name='login', renderer='templates/login.pt')
def user_login(request):
    username = webutils.get(request, 'username', "")
    reason = webutils.get(request, 'reason', None)
    
    return {'layout': site_layout(),
            'username': username,
            'reason':reason,
            'pageTitle': "Login"}


@view_config(route_name='process_login', renderer='templates/information.pt')
def user_login_success(request):
    result = __process_login(request)
    if result == True:
        return {'layout': site_layout(),
                'pageTitle': "Login",
                'header': "Login Successful",
                'message': "You have successfully logged in."}
    else:
        raise HTTPFound(request.application_url+"/login?"+urllib.urlencode(result))

@view_config(route_name='logout', renderer='templates/information.pt')
def user_logout(request):    
    request.response.headers.extend(security.forget(request))
    request.user = None
    
    return {'layout': site_layout(),
            'pageTitle': "Logout",
            'header': "Logout Successful",
            'message': "You have successfully logged out."}

## Internal Functions

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
    username = webutils.post(request, 'username', "")
    pass1    = webutils.post(request, 'pass1', "")
    pass2    = webutils.post(request, 'pass2', "")
    
    name      = webutils.post(request, 'name', "")
    email     = webutils.post(request, 'email', "")
    institute = webutils.post(request, 'institution', "")
    
    resp_dict = {'username': username, 'name':name, 'email':email, 'institution':institute}
    
    if username == "" or pass1 == "" or pass2 == "" or name == "" or email == "" or institute == "":
        resp_dict['reason'] = "Form fields cannot be empty"
        return resp_dict

    if( pass1 != pass2 ):
        resp_dict['reason'] = "Password confirmation does not match"
        return resp_dict

    return True