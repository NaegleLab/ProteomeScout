from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config
import ptmscout.utils.webutils as webutils
from ptmscout.config import strings
from ptmscout.database import experiment, user
from ptmscout.database.user import NoSuchUser
import ptmscout.utils.mail as mail
from ptmscout.database import permissions

@view_config(route_name='privatize_experiment', renderer='ptmscout:templates/accounts/modify_confirm.pt', permission='private')
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
        exp.makePrivate()
        exp.saveExperiment()
        message = strings.privatize_experiment_success_message
        redirect = request.application_url + "/account/experiments"
    
    return {'confirm': confirm,
            'experiment': exp,
            'message': message,
            'redirect': redirect,
            'pageTitle': strings.privatize_experiment_page_title}


@view_config(route_name='publish_experiment', renderer='ptmscout:templates/accounts/modify_confirm.pt', permission='private')
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

@view_config(route_name='invite_experiment', renderer='ptmscout:templates/accounts/modify_confirm.pt', permission='private')
def confirm_invite_user(request):
    confirm = webutils.post(request, 'confirm', "false")
    email = webutils.get(request, 'email', "").strip()
    
    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, request.user)
    
    message = ""
    redirect = None
    
    if email == "":
        message = strings.user_invite_email_required
    elif confirm == "false":
        message = strings.user_invite_confirm % (email)
    else:
        mail.send_automail_message(request, [email], strings.user_invite_email_subject % (request.user.name), strings.user_invite_email_message % (email, request.user.name, exp.name, request.application_url + "/register?email=" + email))
        
        inv = permissions.Invitation(email, eid, request.user.id)
        inv.saveInvitation()
        
        message = strings.user_invited % (email)
        redirect = request.application_url + "/account/experiments"
    
    return {'confirm': confirm,
            'experiment': exp,
            'message': message,
            'redirect': redirect,
            'pageTitle': strings.user_invite_page_title}

@view_config(route_name='share_experiment', renderer='ptmscout:templates/accounts/share.pt', permission='private')
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
                
                mail.send_automail_message(request, [email], strings.user_invite_email_subject % (request.user.name), strings.user_invite_email_message % (new_user.name, request.user.name, exp.name, request.application_url + "/login"))
                
                exp.grantPermission(new_user, 'view')
                exp.saveExperiment()
            except NoSuchUser:
                raise HTTPFound(request.application_url + "/account/experiments/" + str(expid) + "/invite?email=" + email)
        
    return {'pageTitle':strings.share_experiment_page_title,
            'users':users,
            'experiment':exp,
            'reason':reason}

@view_config(route_name='my_experiments', renderer='ptmscout:templates/accounts/my_experiments.pt', permission='private')
def manage_experiments(request):
    users_experiments = [ p.experiment for p in request.user.permissions if p.access_level=='owner' ]
    
    return {'pageTitle':strings.my_experiments_page_title,
            'experiments': users_experiments}
