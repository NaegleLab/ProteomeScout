from pyramid.view import view_config
from ptmscout.layout import site_layout

@view_config(route_name='upload', renderer='templates/information.pt', permission='private')
def user_upload(request):    
    return {'layout': site_layout(),
            'pageTitle': "Upload",
            'header': "Data Upload",
            'message': "Data upload feature is not yet available."}
