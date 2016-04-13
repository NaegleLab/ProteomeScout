from ptmscout.config import strings, settings
from ptmscout.utils import webutils, protein_utils
from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound, HTTPNotFound
from pyramid.response import FileResponse
from ptmscout.database import jobs
from ptmworker import export_tasks
from random import randint
import sys, os, csv, time

@view_config(route_name='batch_download', permission='private')
def download_batch_annotation_file(request):
    fname = "batch.%s.%d.zip" % (request.matchdict['id'], request.user.id)
    fpath = os.path.join(settings.ptmscout_path, settings.annotation_export_file_path, fname)
    if not os.path.exists(fpath):
        raise HTTPNotFound()

    response = FileResponse(fpath, request=request, content_type='application/zip')
    response.content_disposition = 'attachment; filename="%s"' % (fname)
    return response

@view_config(route_name='batch_submit', renderer='ptmscout:templates/info/information.pt', permission='private')
def start_accession_search_job(request):
    return {'pageTitle': strings.protein_batch_search_submitted_page_title,
            'header': strings.protein_batch_search_submitted_page_title,
            'message': strings.protein_batch_search_submitted_message % (request.route_url('my_experiments')) }

def validate_accessions(accessions):
    errors = []
    accessions = accessions.split()
    for acc in accessions:
        if protein_utils.get_accession_type(acc) not in protein_utils.get_valid_accession_types():
            errors.append("Unsupported accession format: %s" % (acc))
    return errors

def create_job_and_submit(request, accessions, user_id):
    accessions = accessions.split()
    batch_id = "%f.%d" % (time.time(), randint(0,10000))

    j = jobs.Job()
    j.status = 'in queue'
    j.stage = 'initializing'
    j.progress = 0
    j.max_progress = 0
    j.status_url = request.route_url('my_experiments')
    j.result_url = request.route_url('batch_download', id=batch_id)
    j.name = "Batch annotate %d proteins" % (len(accessions))
    j.type = 'batch_annotate'
    j.user_id = user_id
    j.save()

    export_tasks.batch_annotate_proteins.apply_async((accessions, batch_id, user_id, j.id))

@view_config(route_name='batch_search', renderer='ptmscout:templates/proteins/batch_search.pt', request_method="POST", permission="private")
def batch_search_view_POST(request):
    terms_of_use_accepted = 'terms_of_use' in request.POST
    accessions = webutils.post(request, 'accessions', "").strip()

    errors = []
    if not terms_of_use_accepted:
        errors = ["You must agree to the ProteomeScout terms of use"]
    elif len(accessions) == 0:
        errors = ["You must enter at least one accession number"]
    else:
        errors = validate_accessions(accessions)

    if len(errors) == 0:
        create_job_and_submit(request, accessions, request.user.id)
        return HTTPFound(request.route_url('batch_submit'))

    return {'pageTitle': strings.protein_batch_search_page_title,
            'errors': errors,
            'accessions': accessions}

@view_config(route_name='batch_search', renderer='ptmscout:templates/proteins/batch_search.pt', request_method="GET", permission="private")
def batch_search_view_GET(request):
    return {'pageTitle': strings.protein_batch_search_page_title,
            'errors': [],
            'accessions':""}
