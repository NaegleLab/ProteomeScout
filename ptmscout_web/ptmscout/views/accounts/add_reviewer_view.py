from pyramid.view import view_config
from ptmscout.database import experiment, user
from ptmscout.utils import webutils, crypto
from ptmscout.config import strings, settings
import urllib

def create_unique_username():
    username = 'reviewer_' + crypto.randomString(4)
    found = False
    while not found:
        try:
            username = 'reviewer_' + crypto.randomString(4)
            user.getUserByUsername(username)
        except user.NoSuchUser:
            found = True
    return username

def create_review_user(exp):
    username = create_unique_username()
    auth_token = crypto.randomString(5)
    
    review_user = user.User(username=username, name="Anonymous Reviewer", email='', institution='', access_level='reviewer')
    review_user.createUser(auth_token)
    review_user.setActive()
    review_user.setExpiration(settings.REVIEWER_ACCOUNT_EXPIRATION_DAYS)
    review_user.saveUser()
    
    exp.grantPermission(review_user, 'view')
    exp.saveExperiment()
    
    return username, auth_token

def get_login_link(exp, request):
    return "%s?%s" % (request.route_url('login'), urllib.urlencode({'redirect':request.route_path('experiment', id=exp.id)}))

@view_config(route_name='review_account', renderer='ptmscout:templates/accounts/modify_confirm.pt')
def create_reviewer_account(request):
    confirm = webutils.post(request, 'confirm', "false") == "true"

    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, request.user)

    message = strings.create_reviewer_account_confirm_message
    header = strings.create_reviewer_account_confirm_header
    redirect = None

    if confirm:
        username, password = create_review_user(exp)
        login_link = get_login_link(exp, request)
        
        message = strings.create_reviewer_account_created_message % (username, password, login_link)
        header = strings.create_reviewer_account_created_header

    return {'confirm': confirm,
            'experiment': exp,
            'message': message,
            'header': header,
            'redirect': redirect,
            'pageTitle': strings.create_reviewer_account_page_title}



