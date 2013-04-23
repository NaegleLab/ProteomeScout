from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config
import ptmscout.utils.webutils as webutils
from ptmscout.config import strings
from ptmscout.database.user import NoSuchUser
import ptmscout.utils.mail as mail
from ptmscout.database import permissions, modifications, experiment, user,\
    upload

@view_config(route_name='privatize_experiment', renderer='ptmscout:templates/accounts/modify_confirm.pt', permission='private')
def privatize_experiment(request):
    confirm = webutils.post(request, 'confirm', "false") == "true"
    
    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, request.user)
    
    message = ""
    header = ""
    redirect = None
    
    if exp.public == 0:
        message = strings.privatize_experiment_already_message
        header = 'Success'
        confirm=True
        redirect = request.application_url + "/account/experiments"
    elif not confirm:
        message = strings.privatize_experiment_confirm_message
        header = 'Confirm'
    else:
        exp.makePrivate()
        exp.saveExperiment()
        message = strings.privatize_experiment_success_message
        header = 'Success'
        redirect = request.application_url + "/account/experiments"
    
    return {'confirm': confirm,
            'experiment': exp,
            'message': message,
            'header': header,
            'redirect': redirect,
            'pageTitle': strings.privatize_experiment_page_title}


@view_config(route_name='publish_experiment', renderer='ptmscout:templates/accounts/modify_confirm.pt', permission='private')
def publish_experiment(request):
    confirm = webutils.post(request, 'confirm', "false") == "true"
    
    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, request.user)
    
    message = ""
    header = ""
    redirect = None
    
    if exp.public == 1:
        message = strings.publish_experiment_already_message
        header = "Success"
        confirm=True
        redirect = request.application_url + "/account/experiments"
    elif not confirm:
        message = strings.publish_experiment_confirm_message
        header = "Confirm"
    else:
        exp.makePublic()
        exp.saveExperiment()
        message = strings.publish_experiment_success_message
        header = "Success"
        redirect = request.application_url + "/account/experiments"
    
    return {'confirm': confirm,
            'experiment': exp,
            'header': header,
            'message': message,
            'redirect': redirect,
            'pageTitle': strings.publish_experiment_page_title}

@view_config(route_name='invite_experiment', renderer='ptmscout:templates/accounts/modify_confirm.pt', permission='private')
def confirm_invite_user(request):
    confirm = webutils.post(request, 'confirm', "false") == "true"
    email = webutils.get(request, 'email', "").strip()
    
    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, request.user)
    
    message = ""
    header = ""
    redirect = request.application_url + "/account/experiments"
    
    email_ok, domain_ok = mail.email_is_valid(email)
    
    if email == "":
        message = strings.user_invite_email_required
        header = "Failure"
        confirm = True
    elif not email_ok:
        message = strings.failure_reason_email_not_valid
        header = "Failure"
        confirm = True
    elif not domain_ok:
        message = strings.failure_reason_email_not_allowed
        header = "Failure"
        confirm = True
    elif not confirm:
        message = strings.user_invite_confirm % (email)
        header = "Confirm"
        redirect = None
    else:
        mail.send_automail_message(request, [email], strings.user_invite_email_subject % (request.user.name), strings.user_invite_email_message % (email, request.user.name, exp.name, request.application_url + "/register?email=" + email))
        
        inv = permissions.Invitation(email, eid, request.user.id)
        inv.saveInvitation()
        
        message = strings.user_invited % (email)
        header = "Success"
    
    return {'confirm': confirm,
            'experiment': exp,
            'header': header,
            'message': message,
            'redirect': redirect,
            'pageTitle': strings.user_invite_page_title}


def add_user_permission(request, users, exp):

    email = webutils.post(request, 'email', "").strip()
    if email == "":
        return strings.failure_reason_form_fields_cannot_be_empty
    else:
        try:
            new_user = user.getUserByEmail(email)
            users.append(new_user)
            
            mail.send_automail_message(request, [email], strings.user_invite_email_subject % (request.user.name), strings.user_invite_email_message % (new_user.name, request.user.name, exp.name, request.application_url + "/login"))
            
            exp.grantPermission(new_user, 'view')
            exp.saveExperiment()
        except NoSuchUser:
            raise HTTPFound(request.application_url + "/account/experiments/" + str(exp.id) + "/invite?email=" + email)
    return None

def remove_user_permission(request, users, exp):
    uid = webutils.post(request, 'uid', None)
    try:
        uid = int(uid)
        removed_user = None
        for u in users:
            if u.id == uid:
                removed_user = u
                break
        users.remove(removed_user)
        exp.revokePermission(removed_user)
        exp.saveExperiment()
    except:
        pass

@view_config(route_name='share_experiment', renderer='ptmscout:templates/accounts/share.pt', permission='private')
def manage_experiment_permissions(request):
    reason = None
    submitted = webutils.post(request, 'submitted', "0")
    mode = webutils.post(request, 'mode', "").strip()
    
    expid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(expid, request.user)
    
    users = [ p.user for p in exp.permissions if p.user.id != request.user.id ]
    
    if(submitted == "1"):
        if mode == 'add':
            reason = add_user_permission(request, users, exp)
        elif mode == 'remove':
            reason = remove_user_permission(request, users, exp)

    return {'pageTitle':strings.share_experiment_page_title,
            'users':users,
            'experiment':exp,
            'reason':reason}

def get_sessions(user_experiments):
    session_map = {}
    
    for exp in user_experiments:
        session = upload.getMostRecentSession(exp.id)
        session_map[exp.id] = session.id
    
    return session_map

def get_experiment_breakdown( users_experiments ):
    in_process = [ exp for exp in users_experiments if exp.status == 'configuration' ]
    available = [ exp for exp in users_experiments if exp.status != 'configuration' ]
    error_state = [ exp for exp in users_experiments if exp.status == 'error' ]
    sessions = get_sessions(in_process + error_state)
    
    return {
            'available': available,
            'in_process': in_process,
            'sessions': sessions,
           }

@view_config(route_name='my_experiments', renderer='ptmscout:templates/accounts/my_experiments.pt', permission='private')
def manage_experiments(request):
    users_experiments = request.user.myExperiments()
    users_datasets = request.user.myDatasets()

    users_active_jobs = [ job for job in request.user.jobs if job.is_active() or not job.is_old() ]
    
    
    return { 'pageTitle':strings.my_experiments_page_title,
             'experiments':get_experiment_breakdown( users_experiments ),
             'datasets':get_experiment_breakdown( users_datasets ),
             'jobs': users_active_jobs }
