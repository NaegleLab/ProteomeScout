from ptmscout.config import strings
from ptmscout.database import user
from ptmscout.database.user import NoSuchUser
from ptmscout.utils import webutils
from pyramid.view import view_config

@view_config(route_name='activate_account', renderer='ptmscout:templates/info/information.pt')
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
    
    url_redirect = request.application_url+"/login"
    return {'pageTitle': strings.account_activation_page_title,
            'header': strings.account_activation_success_header,
            'message': strings.account_activation_success_message % (url_redirect),
            'redirect': url_redirect}


## Internal Functions

def __failed_account_activation(request):
    url_redirect = request.application_url+"/register"
    return {'pageTitle':strings.account_activation_page_title,
            'header':strings.account_activation_failed_header,
            'message':strings.account_activation_failed_message % (url_redirect),
            'redirect': url_redirect}
