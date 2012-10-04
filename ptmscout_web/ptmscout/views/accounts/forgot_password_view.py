from ptmscout.config import strings
from ptmscout.database import user
from ptmscout.database.user import NoSuchUser
from ptmscout.utils import webutils, crypto, mail
from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config
import urllib

@view_config(route_name='forgot_password', renderer='ptmscout:templates/accounts/forgot_password.pt')
def forgot_password(request):
    email = webutils.get(request, 'email', "")
    reason = webutils.get(request, 'reason', None)
    
    return {'pageTitle': strings.forgotten_password_page_title,
            'email': email,
            'reason': reason}

@view_config(route_name='process_forgot_password', renderer='ptmscout:templates/info/information.pt')
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
            'message': strings.forgotten_password_success_message,
            'redirect': request.application_url+  "/login"}
