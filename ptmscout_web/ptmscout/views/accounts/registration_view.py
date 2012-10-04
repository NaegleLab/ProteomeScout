from ptmscout.config import settings, strings
from ptmscout.database import user
from ptmscout.database.user import NoSuchUser
from ptmscout.utils import webutils, mail
from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config
import re
import urllib




@view_config(route_name='register',renderer='ptmscout:templates/accounts/user_registration.pt')
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
    
@view_config(route_name='process_registration', renderer='ptmscout:templates/info/information.pt')
def user_registration_success(request):
    result = __process_registration(request)
    if result == True:
        return {
            'pageTitle': strings.user_registration_page_title,
            'header': strings.user_registration_success_header,
            'message': strings.user_registration_success_message,
            'redirect': request.application_url + "/login"}
    else:
        raise HTTPFound(request.application_url+"/register?"+urllib.urlencode(result))


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

    if len(pass1) < settings.MINIMUM_PASSWORD_LENGTH:
        resp_dict['reason'] = strings.failure_reason_password_too_short % settings.MINIMUM_PASSWORD_LENGTH
        return resp_dict

    if( pass1 != pass2 ):
        resp_dict['reason'] = strings.failure_reason_new_passwords_not_matching
        return resp_dict

    ptm_user = user.User(username, name, email, institute)
    ptm_user.createUser(pass1)
    
    ptm_user.processInvitations()
    
    ptm_user.saveUser()
    
    mail.send_automail_message(request, 
    [email], 
    strings.user_registration_email_subject, 
    strings.user_registration_email_message % (ptm_user.name, request.application_url, ptm_user.username, ptm_user.activation_token))
    
    return True