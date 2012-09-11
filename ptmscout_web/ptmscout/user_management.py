from pyramid.httpexceptions import HTTPFound, HTTPForbidden
from pyramid.view import view_config
import urllib
import utils.webutils as webutils
from utils import crypto
import strings
from ptmscout.database import experiment, user
from ptmscout.database.user import NoSuchUser

@view_config(route_name='privatize_experiment', renderer='templates/modify_confirm.pt', permission='private')
def privatize_experiment(request):
    confirm = webutils.post(request, 'confirm', "false")
    
    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, request.user)
    
    message = ""
    redirect = None
    
    if exp.public == 0:
        message = strings.privatize_experiment_already_message
        confirm='true'
        redirect = request.application_url + "/account/experiments"
    elif confirm == "false":
        message = strings.privatize_experiment_confirm_message
    else:
        exp.makePublic()
        exp.saveExperiment()
        message = strings.privatize_experiment_success_message
        redirect = request.application_url + "/account/experiments"
    
    return {'confirm': confirm,
            'experiment': exp,
            'message': message,
            'redirect': redirect,
            'pageTitle': strings.privatize_experiment_page_title}


@view_config(route_name='publish_experiment', renderer='templates/modify_confirm.pt', permission='private')
def publish_experiment(request):
    confirm = webutils.post(request, 'confirm', "false")
    
    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, request.user)
    
    message = ""
    redirect = None
    
    if exp.public == 1:
        message = strings.publish_experiment_already_message
        confirm='true'
        redirect = request.application_url + "/account/experiments"
    elif confirm == "false":
        message = strings.publish_experiment_confirm_message
    else:
        exp.makePublic()
        exp.saveExperiment()
        message = strings.publish_experiment_success_message
        redirect = request.application_url + "/account/experiments"
    
    return {'confirm': confirm,
            'experiment': exp,
            'message': message,
            'redirect': redirect,
            'pageTitle': strings.publish_experiment_page_title}

@view_config(route_name='share_experiment', renderer='templates/share.pt', permission='private')
def manage_experiment_permissions(request):
    reason = None
    submitted = webutils.post(request, 'submitted', "0")
    email = webutils.post(request, 'email', "").strip()
    
    expid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(expid, request.user)
    
    users = [ p.user for p in exp.permissions if p.user.id != request.user.id ]
    
    if(submitted == "1"):
        if email == "":
            reason = strings.failure_reason_form_fields_cannot_be_empty
        else:
            try:
                new_user = user.getUserByEmail(email)
                users.append(new_user)
                exp.grantPermission(new_user, 'view')
                exp.saveExperiment()
            except NoSuchUser:
                reason = strings.failure_reason_email_address_not_on_record
        
    return {'pageTitle':strings.share_experiment_page_title,
            'users':users,
            'experiment':exp,
            'reason':reason}

@view_config(route_name='my_experiments', renderer='templates/my_experiments.pt', permission='private')
def manage_experiments(request):
    users_experiments = [ p.experiment for p in request.user.permissions if p.access_level=='owner' ]
    
    return {'pageTitle':strings.my_experiments_page_title,
            'experiments': users_experiments}

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
            'header': strings.success_header,
            'redirect': request.application_url + "/account"}
