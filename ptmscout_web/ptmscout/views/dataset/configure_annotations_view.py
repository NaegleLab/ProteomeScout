from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.utils import uploadutils, webutils, decorators
from ptmscout.database import upload, experiment, modifications
import re
from pyramid.httpexceptions import HTTPFound

def check_is_type(val, tp):
    try:
        tp(val)
        return True
    except:
        return False

def verify_dataset(session, data_rows, valid_MS_ids):
    errors = []
    
    for i, row in enumerate(data_rows):
        for sc in session.columns:
            pe = None
            if sc.type != 'none' and sc.column_number >= len(row):
                pe = uploadutils.ParseError(i+1, sc.column_number, strings.experiment_upload_warning_missing_column)
                
            elif sc.type == 'MS_id' and not check_is_type(row[sc.column_number], int):
                pe = uploadutils.ParseError(i+1, sc.column_number, strings.experiment_upload_warning_columns_values_should_be % ("integer"))

            elif sc.type == 'MS_id' and not int(row[sc.column_number]) in valid_MS_ids:
                pe = uploadutils.ParseError(i+1, sc.column_number, strings.experiment_upload_error_ms_id_not_found % int(row[sc.column_number]))

            elif sc.type == 'numeric' and not check_is_type(row[sc.column_number], float):
                pe = uploadutils.ParseError(i+1, sc.column_number, strings.experiment_upload_warning_columns_values_should_be % ("numeric"))
            
            if pe != None:
                errors.append(pe)
    
    return [ e.message for e in errors ]

def parse_column_assignments(request, session, headers):
    column_settings = {}
    
    for field in request.POST:
        m = re.match('column_([0-9]+)_(type|label)', field)
        if m:
            col = int(m.group(1))
            col_def = column_settings.get(col, {'type':'','label':''})
            col_def[m.group(2)] = request.POST[field].strip()
            
            column_settings[col] = col_def
    
    errors = []
    
    found_MScolumn = False
    found_annotation = False
    label_required = set(['numeric', 'nominative', 'cluster'])
    
    session.columns = []
    
    for c in column_settings:
        col_type = column_settings[c]['type']
        col_label = column_settings[c]['label']
        
        sc = upload.SessionColumn()
        sc.column_number = c
        sc.type = col_type
        sc.label = col_label
        
        if col_type == 'MS_id':
            if found_MScolumn:
                errors.append(strings.experiment_upload_error_limit_one_column_of_type % ("MS_id"))
            found_MScolumn = True
        elif col_type in label_required:
            if col_label == '':
                errors.append(strings.experiment_upload_error_data_column_empty_label % (c) )
            found_annotation = True
        
        session.columns.append(sc)
        
    if not found_MScolumn:
        errors.append(strings.experiment_upload_error_no_column_assignment % ('MS_id'))
        
    if not found_annotation:
        errors.append(strings.experiment_upload_error_no_annotations)       
    
    session.stage = 'confirm'
    session.save()
    
    return errors

def configure_annotations_POST(request, experiment, session):
    force = webutils.post(request, 'override', "false") != "false"
    measurements = modifications.getMeasuredPeptidesByExperiment(experiment.id, request.user)
    valid_MS_ids = set([ms.id for ms in measurements])

    headers, data_rows = uploadutils.load_header_and_data_rows(session.data_file, N=100)
    errors = parse_column_assignments(request, session, headers)
    if len(errors) > 0:
        result = configure_annotations_GET(request, experiment, session)
        result['error'] = errors
        return result
        
    errors = verify_dataset(session, data_rows, valid_MS_ids)
    
    if len(errors) > 0 and not force:
        result = configure_annotations_GET(request, experiment, session)
        result['error'] = errors
        result['allowoverride'] = True
        return result
    
    return HTTPFound(request.route_url('confirm_annotations', id=experiment.id, sid=session.id))

@view_config(route_name='configure_annotations', request_method='POST', renderer='ptmscout:/templates/upload/upload_config.pt', permission='private')
@decorators.get_experiment('id',types=set(['experiment']))    
@decorators.get_session('sid','annotations')
def configure_annotations_POST_view(context, request, experiment, session):
    return configure_annotations_POST(request, experiment, session)
    
    
def get_columns_from_header(h):
    col = {'type':'','label':''}
    
    if h == "MS_id":
        col['type'] = 'MS_id'
    
    m = re.match('^([a-zA-Z]+):(.+)$', h)
    if(m and m.group(1) == 'cluster'):
        col['type'] = 'cluster'
        col['label'] = m.group(2)
    elif(m and m.group(1) == 'numeric'):
        col['type'] = 'numeric'
        col['label'] = m.group(2)
    elif(m and m.group(1) == 'nominative'):
        col['type'] = 'nominative'
        col['label'] = m.group(2)
    
    return col



def configure_annotations_GET(request, experiment, session):
    headers, data_rows = uploadutils.load_header_and_data_rows(session.data_file, N=100, truncate=30)
    
    data_definitions = { 'columns':dict([ (i, get_columns_from_header(h)) for i, h in enumerate(headers) ] ) }
    if len(session.columns) == len(headers):
        session_defs = uploadutils.assign_columns_from_session(session)
        session_defs.update(data_definitions)
        data_definitions = session_defs
        
    return {
            'session_id': session.id,
            'pageTitle':strings.upload_configure_annotations_page_title,
            'error':[],
            'headers':headers,
            'data_rows':data_rows,
            'allowoverride':False,
            'data_definitions': data_definitions,
            'column_values': ['none','MS_id','cluster','numeric','nominative']
            }

@view_config(route_name='configure_annotations', request_method='GET', renderer='ptmscout:/templates/upload/upload_config.pt', permission='private')
@decorators.get_experiment('id',types=set(['experiment']))    
@decorators.get_session('sid','annotations')
def configure_annotations_GET_view(context, request, experiment, session):
    return configure_annotations_GET(request, experiment, session)
