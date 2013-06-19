from pyramid.view import view_config
from ptmscout.utils import decorators
from pyramid.response import FileResponse
from pyramid.httpexceptions import HTTPForbidden
from ptmscout.utils import webutils
from ptmscout.database import experiment, jobs
from ptmscout.utils import downloadutils
from ptmscout.config import strings, settings
from ptmworker import export_tasks
import time
import os

@view_config(route_name='experiment_download', renderer='tsv', permission='private')
def download_experiment(request):
    get_errors = bool( webutils.get(request, 'errors', False) )
    exp_id = int(request.matchdict['id'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    if exp not in request.user.myExperiments():
        raise HTTPForbidden()

    header, rows = downloadutils.annotate_experiment_with_errors(exp, get_errors)
    
    exp_filename = 'experiment.%d.tsv' % (exp_id)
    if get_errors:
        exp_filename = 'experiment.%d.errors.tsv' % (exp_id)
   
    request.response.content_type = 'text/tab-separated-values'
    request.response.content_disposition = 'attachment; filename="%s"' % (exp_filename)
    return { 'header': header, 'data': rows }


@view_config(route_name='export_download', permission='private')
@decorators.get_experiment('id')
def download_export_file(context, request, exp):
    fname = "experiment.%d.%d.%s.tsv" % (exp.id, request.user.id, request.matchdict['export_id'])
    
    fpath = os.path.join(settings.ptmscout_path, settings.annotation_export_file_path, fname)
    if not os.path.exists(fpath):
        raise HTTPNotFound()

    response = FileResponse(fpath, request)
    request.response.content_type = 'text/tab-separated-values'
    request.response.content_disposition = 'attachment; filename="%s"' % (fname)
    return response


def create_annotate_export_job(request, export_id, exp_id, user_id):
    export_job = jobs.Job()
    export_job.name = "Export Experiment %d" % (exp_id)
    export_job.user_id = user_id
    export_job.result_url = request.route_url('export_download', id=exp_id, export_id=export_id)
    export_job.status_url = request.route_url('my_experiments')
    export_job.status = 'in queue'
    export_job.type = 'export_experiment'
    
    export_job.save()

    return export_job.id
    
@view_config(route_name='experiment_export', renderer="ptmscout:templates/info/information.pt", permission='private')
@decorators.get_experiment('id', types=set(['experiment','dataset']))
def export_experiment(context, request, exp):
    annotate = webutils.get(request, 'annotate', 'no') == 'yes'
    user_id = request.user.id
    exp_id = exp.id
    export_id = int(time.time())

    job_id = create_annotate_export_job(request, export_id, exp_id, user_id)
    export_tasks.run_experiment_export_job.apply_async( (annotate, export_id, exp_id, user_id, job_id) )

    return { 'pageTitle':strings.experiment_export_started_page_title,
             'header':strings.experiment_export_started_page_title,
             'message': strings.experiment_export_started_message % (request.route_url('my_experiments')) }
