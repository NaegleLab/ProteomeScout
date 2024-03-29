from ptmscout.config import strings
from pyramid.view import view_config
from ptmscout.views.upload import upload_configure
from ptmscout.utils import decorators, wizard

def create_nav_wizard(request, session):
    navigation = wizard.WizardNavigation(request)

    navigation.add_page('dataset_configure', "Configure Dataset", True, id=session.id)
    navigation.add_page('dataset_confirm', "Confirm Upload", False, id=session.id)
    navigation.set_page('dataset_configure')

    return navigation


def upload_configure_dataset(request, session):
    rval = upload_configure.upload_config_handler( request, session, strings.dataset_upload_configure_page_title, create_nav_wizard(request, session), mod_required=False, nextStage='confirm')

    if isinstance(rval, dict):
        rval['instructions'] = strings.dataset_upload_configure_instructions

    return rval

@view_config(route_name='dataset_configure', renderer='ptmscout:/templates/upload/upload_config.pt')
@decorators.get_session('id', 'dataset')
def upload_configure_dataset_view(context, request, session):
    return upload_configure_dataset(request, session)
