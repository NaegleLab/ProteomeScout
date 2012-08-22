from layout import site_layout
from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config
import urllib
import webutils



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


## Internal Functions

def __process_login(request):
    username = webutils.post(request, 'username', "")
    password = webutils.post(request, 'password', "")
    
    resp_dict = {'username': username}
    
    if username == "" or password == "":
        resp_dict['reason'] = "All fields are required"
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