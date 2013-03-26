from pyramid.view import view_config

@view_config(route_name='review_account', renderer='ptmscout:templates/accounts/modify_confirm.pt')
def create_reviewer_account(request):
    confirm = webutils.post(request, 'confirm', "false") == "true"

    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, request.user)

    message = strings.create_reviewer_account_confirm_message
    header = strings.create_reviewer_account_confirm_header
    redirect = None

    if confirm:
        message = strings.create_reviewer_account_created_message
        header = strings.create_reviewer_account_created_header

    return {'confirm': confirm,
            'experiment': exp,
            'message': message,
            'header': header,
            'redirect': redirect,
            'pageTitle': strings.create_reviewer_account_page_title}



