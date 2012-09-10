from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config
import urllib
import utils.webutils as webutils
from utils import crypto
import strings

@view_config(route_name='account_management', renderer='templates/account.pt', permission='private')
def manage_account(request):
    reason = webutils.get(request, 'reason', None)
    return {'username': request.user.username,
            'fullname': request.user.name,
            'email': request.user.email,
            'institution': request.user.institution,
            'pageTitle': strings.account_management_page_title,
            'reason': reason}

@view_config(route_name='change_password', permission='private')
def change_password(request):
    oldpass = webutils.post(request, 'old_pass', "")
    newpass1 = webutils.post(request, 'new_pass1', "")
    newpass2 = webutils.post(request, 'new_pass2', "")
    
    if oldpass == "" or newpass1 == "" or newpass2 == "":
        raise HTTPFound(request.application_url + "/account?" + urllib.urlencode({'reason': strings.failure_reason_form_fields_cannot_be_empty}))
    
    if newpass1 != newpass2:
        raise HTTPFound(request.application_url + "/account?" + urllib.urlencode({'reason': strings.failure_reason_new_passwords_not_matching}))
    
    _, old_salted_pass = crypto.saltedPassword(oldpass, request.user.salt)
    
    if old_salted_pass != request.user.salted_password:
        raise HTTPFound(request.application_url + "/account?" + urllib.urlencode({'reason': strings.failure_reason_incorrect_password}))
    
    _, new_salted_pass = crypto.saltedPassword(newpass1, request.user.salt)
    request.user.salted_password = new_salted_pass
    request.user.saveUser()
    
    return HTTPFound(request.application_url + "/change_password_success")
    
@view_config(route_name='change_password_success', renderer='templates/information.pt', permission='private')
def change_password_success(request):
    return {'pageTitle': strings.change_password_page_title,
            'message': strings.change_password_success_message,
            'header': strings.success_header}
