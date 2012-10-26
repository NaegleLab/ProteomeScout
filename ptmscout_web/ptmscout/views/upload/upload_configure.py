from ptmscout.config import strings
from pyramid.view import view_config


@view_config(route_name='upload_config', renderer='ptmscout:/templates/upload/upload_config.pt')
def upload_config(request):
    return {'pageTitle': strings.experiment_upload_configure_page_title,
            'instruction': strings.experiment_upload_configure_message}