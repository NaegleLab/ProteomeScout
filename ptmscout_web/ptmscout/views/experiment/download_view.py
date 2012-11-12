from pyramid.view import view_config
from ptmscout.utils import uploadutils, webutils
import sys
from ptmscout.database import experiment
from pyramid.httpexceptions import HTTPForbidden
from ptmscout.config import strings

def insert_errors(errors, rows):
    for row in rows:
        row.insert(0, [])
    
    for error in errors:
        rows[error.line-1][0].append(error.message)
        
    for row in rows:
        row[0] = ", ".join(row[0])

@view_config(route_name='experiment_download', renderer='tsv', permission='private')
def download_experiment(request):
    get_errors = bool( webutils.get(request, 'errors', False) )
    exp_id = int(request.matchdict['id'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    if exp not in request.user.myExperiments():
        raise HTTPForbidden()
    
    header, rows = uploadutils.load_header_and_data_rows(exp.dataset, sys.maxint)
    
    if get_errors:
        header.insert(0, strings.experiment_upload_error_reasons_column_title)
        insert_errors(exp.errors, rows)
        
        rows = [r for r in rows if r[0] != ""]
    
    return { 'header': header, 'data': rows }