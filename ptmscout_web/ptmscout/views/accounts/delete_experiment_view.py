from ptmscout.utils import decorators, webutils
from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.database import modifications


@view_config(route_name='delete_experiment', renderer='ptmscout:/templates/accounts/modify_confirm.pt', permission="private", request_method="POST")
@decorators.get_experiment('id', types=set(['dataset']), owner_required=True)
def view_function_POST(context, request, exp):
    modifications.deleteExperimentData(exp.id)
    exp.delete()
    
    message = strings.delete_experiment_success_message
    header = 'Success'
    redirect = request.route_url('my_experiments')
    
    return {'confirm': True,
            'experiment': exp,
            'message': message,
            'header': header,
            'redirect': redirect,
            'pageTitle': strings.delete_experiment_page_title}
    
@view_config(route_name='delete_experiment', renderer='ptmscout:/templates/accounts/modify_confirm.pt', permission="private", request_method="GET")
@decorators.get_experiment('id', types=set(['dataset']), owner_required=True)
def view_function_GET(context, request, exp):
    message = strings.delete_experiment_confirm_message
    header = 'Confirm'

    return {'confirm': False,
            'experiment': exp,
            'message': message,
            'header': header,
            'redirect': None,
            'pageTitle': strings.delete_experiment_page_title}