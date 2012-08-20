from database.user import PTMUser
from layout import site_layout
from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config
import urllib

def process_registration(request):
    try:
        username = request.POST['username']
        pass1 = request.POST['pass1']
        pass2 = request.POST['pass2']
        
        name = request.POST['name']
        email = request.POST['email']
        institute = request.POST['institution']
        
        resp_dict = {"username": username, "name":name, "email":email, "institution":institute}
        
        if username == "" or pass1 == "" or pass2 == "" or name == "" or email == "" or institute == "":
            resp_dict["reason"] = "Form fields cannot be empty"
            return resp_dict

        if( pass1 != pass2 ):
            resp_dict["reason"] = "Password confirmation does not match"
            return resp_dict
        
    except KeyError:
        return {"reason": "Form fields cannot be empty"}
        
    return True

@view_config(route_name='register',renderer='templates/user_registration.pt')
def user_registration_view(request):
    username = ""
    email = ""
    name = ""
    institution = ""
    reason = ""
    
    if "username" in request.GET:
        username = request.GET["username"]
    if "email" in request.GET:
        email = request.GET["email"]
    if "name" in request.GET:
        name = request.GET["name"]
    if "institution" in request.GET:
        institution = request.GET["institution"]
    
    if "reason" in request.GET:
        reason = request.GET["reason"]
    
    return {'layout': site_layout(),
            "username":username,
            "email":email,
            "name":name,
            "institution":institution,
            "reason":reason,
            'pageTitle': "User Registration"}

@view_config(route_name='process_registration', renderer='templates/information.pt')
def user_registration_success(request):
    result = process_registration(request)
    if result == True:
        return {'layout': site_layout(),
            'pageTitle': "User Registration",
            'header': "Registration Successful",
            'message': "A confirmation e-mail has been sent to the specified e-mail address. Please check your e-mail to complete your registration."}
    else:
        raise HTTPFound(request.application_url+"/register?"+urllib.urlencode(result))
