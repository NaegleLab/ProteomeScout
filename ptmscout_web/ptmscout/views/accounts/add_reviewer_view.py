from pyramid.view import view_config
from ptmscout.database import experiment, user
from ptmscout.utils import webutils, crypto
from ptmscout.config import strings



def create_unique_username():
    username = 'reviewer_' + crypto.randomString(5)
    found = False
    while not found:
        try:
            username = 'reviewer_' + crypto.randomString(5)
            user.getUserByUsername(username)
        except user.NoSuchUser:
            found = True
    return username

def create_review_user(exp):
    username = create_unique_username()
    auth_token = crypto.randomString(10)
    
    review_user = user.User(username=username, name="Anonymous Reviewer", email='', institution='')
    review_user.createUser(auth_token)
    review_user.setActive()
    review_user.saveUser()
    
    exp.grantPermission(review_user, 'view')
    exp.saveExperiment()

@view_config(route_name='review_account', renderer='ptmscout:templates/accounts/modify_confirm.pt')
def create_reviewer_account(request):
    confirm = webutils.post(request, 'confirm', "false") == "true"

    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, request.user)

    message = strings.create_reviewer_account_confirm_message
    header = strings.create_reviewer_account_confirm_header
    redirect = None

    if confirm:
        login_link = create_review_user(exp)
        
        message = strings.create_reviewer_account_created_message % (login_link)
        header = strings.create_reviewer_account_created_header

    return {'confirm': confirm,
            'experiment': exp,
            'message': message,
            'header': header,
            'redirect': redirect,
            'pageTitle': strings.create_reviewer_account_page_title}



