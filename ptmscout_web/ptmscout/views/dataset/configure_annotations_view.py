from pyramid.view import view_config
from ptmscout.config import strings



@view_config(route_name='configure_annotations', request_method='POST', renderer='ptmscout:/templates/upload/upload_config.pt', permission='private')
def configure_annotations_POST(request):
    return configure_annotations_GET(request)

@view_config(route_name='configure_annotations', request_method='GET', renderer='ptmscout:/templates/upload/upload_config.pt', permission='private')
def configure_annotations_GET(request):
    experiment_id = int(request.matchdict['id'])
    session_id = int(request.matchdict['sid'])
    
    data_definitions = {'columns':[]}
    headers = []
    data_rows = []
    
    return {
            'session_id': session_id,
            'pageTitle':strings.upload_configure_annotations_page_title,
            'error':[],
            'headers':headers,
            'data_rows':data_rows,
            'allowoverride':False,
            'data_definitions': data_definitions
            }