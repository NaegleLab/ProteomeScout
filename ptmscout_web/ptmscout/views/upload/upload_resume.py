from ptmscout.database import upload
from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound


@view_config(route_name='upload_resume', permission='private')
def resume_upload_session(request):
    session_id = int(request.matchdict['id'])
    
    session = upload.getSessionById(session_id, request.user)
    
    base_url = request.application_url + "/upload/%d/%s"
    
    if session.stage == 'complete':
        return HTTPFound(base_url % (session_id, 'confirm'))
    elif session.stage == 'condition':
        return HTTPFound(base_url % (session_id, 'conditions'))
    else:
        return HTTPFound(base_url % (session_id, session.stage))