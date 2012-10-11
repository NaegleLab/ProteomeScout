from pyramid.view import view_config
from ptmscout.views.layout import site_layout
from ptmscout.config import strings
import base64
import json

def experiment_to_dict(exp):
    expd = {}
    
    for key in exp.__dict__:
        if key[0] != "_":
            expd[key] = exp.__dict__[key]
    
    expd['URL'] = exp.getUrl()
    
    return expd

@view_config(route_name='upload', renderer='ptmscout:templates/accounts/upload.pt', permission='private')
def user_upload(request):    
    users_experiments = [ p.experiment for p in request.user.permissions if p.access_level=='owner' ]
    dict_exps = [ experiment_to_dict(exp) for exp in users_experiments ]
    
    return {'pageTitle': strings.upload_page_title,
            'user_experiments': dict_exps,
            'json_user_data': base64.b64encode(json.dumps(dict_exps))}
