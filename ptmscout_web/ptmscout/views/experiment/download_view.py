from pyramid.view import view_config
from ptmscout.utils import webutils
from ptmscout.database import experiment
from pyramid.httpexceptions import HTTPForbidden
from ptmscout.utils import downloadutils

@view_config(route_name='experiment_download', renderer='tsv', permission='private')
def download_experiment(request):
    get_errors = bool( webutils.get(request, 'errors', False) )
    exp_id = int(request.matchdict['id'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    if exp not in request.user.myExperiments():
        raise HTTPForbidden()

    header, rows = downloadutils.annotate_experiment_file(exp, get_errors)
   
    request.response.content_type = 'text/tab-separated-values'
    request.response.content_disposition = 'attachment; filename="experiment.%d.annotated.tsv"' % (exp_id)
    return { 'header': header, 'data': rows }
