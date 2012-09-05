from pyramid.view import view_config
from ptmscout.layout import site_layout
import strings

@view_config(route_name='upload', renderer='templates/information.pt', permission='private')
def user_upload(request):    
    return {'layout': site_layout(),
            'pageTitle': strings.upload_page_title,
            'header': strings.upload_page_header,
            'message': "Data upload feature is not yet available."}
