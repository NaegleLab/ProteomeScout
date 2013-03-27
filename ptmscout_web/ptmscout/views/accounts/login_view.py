from ptmscout.config import strings
from ptmscout.database import user
from ptmscout.database.user import NoSuchUser
from ptmscout.utils import crypto, webutils
from pyramid import security
from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config
import urllib



@view_config(route_name='login', renderer='ptmscout:templates/accounts/login.pt')
def user_login(request):
    username = webutils.get(request, 'username', "")
    reason = webutils.get(request, 'reason', None)
    redirect = webutils.get(request, 'redirect', None)
    
    return {
            'username': username,
            'reason':reason,
            'pageTitle': strings.login_page_title,
            'origin': redirect
            }


@view_config(route_name='process_login', renderer='ptmscout:templates/info/information.pt')
def user_login_success(request):
    redirect = webutils.post(request, 'origin', None)
    
    result = __process_login(request)
    if result == True:
        if redirect == None:
            redirect = request.application_url
        else:
            redirect = request.application_url + redirect
        
        return {
                'pageTitle': strings.login_page_title,
                'header': strings.login_page_success_header,
                'message': strings.login_page_success_message,
                'redirect': redirect
                }
    else:
        if redirect != None:
            result['redirect'] = redirect
        raise HTTPFound(request.application_url+"/login?"+urllib.urlencode(result))

@view_config(route_name='logout', renderer='ptmscout:templates/info/information.pt')
def user_logout(request):    
    request.response.headers.extend(security.forget(request))
    request.user = None
    
    return {
            'pageTitle': strings.logout_page_title,
            'header': strings.logout_page_header,
            'message': strings.logout_page_message,
            'redirect': request.application_url}
    
    
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
        request.username = username
        request.response.headers.extend(security.remember(request, ptm_user.username))
    except NoSuchUser:
        resp_dict['reason'] = strings.failure_reason_incorrect_credentials
        return resp_dict
    
    return True
