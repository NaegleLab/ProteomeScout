from ptmscout.config import strings
from pyramid.view import view_config
from ptmscout.views.upload import upload_configure


@view_config(route_name='dataset_configure', renderer='ptmscout:/templates/upload/upload_config.pt')
def upload_configure_dataset(request):
    session_id = int(request.matchdict['id'])
    return upload_configure.upload_config_handler( request, strings.dataset_upload_configure_page_title, request.route_url('dataset_confirm', id=session_id), mod_required=False, nextStage='confirm' )
