from ptmscout.config import strings
from pyramid.view import view_config
from ptmscout.views.upload import upload_configure
from ptmscout.utils import decorators

def upload_configure_dataset(request, session):
    return upload_configure.upload_config_handler( request, session, strings.dataset_upload_configure_page_title, request.route_url('dataset_confirm', id=session.id), mod_required=False, nextStage='confirm' )

@view_config(route_name='dataset_configure', renderer='ptmscout:/templates/upload/upload_config.pt')
@decorators.get_session('id', 'dataset')
def upload_configure_dataset_view(context, request, session):
    return upload_configure_dataset(request, session)