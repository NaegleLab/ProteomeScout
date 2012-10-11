from pyramid.view import view_config
from ptmscout.views.layout import site_layout
from ptmscout.config import strings

@view_config(route_name='upload', renderer='ptmscout:templates/accounts/upload.pt', permission='private')
def user_upload(request):    
    users_experiments = [ p.experiment for p in request.user.permissions if p.access_level=='owner' ]
    
    return {'pageTitle': strings.upload_page_title,
            'user_experiments': [(e.id, e.name) for e in users_experiments]}
