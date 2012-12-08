from ptmscout.config import strings
from pyramid.view import view_config
from ptmscout.database import upload
from ptmscout.utils import webutils
from pyramid.httpexceptions import HTTPFound
from ptmscout.utils import uploadutils
import re


def parse_user_input(session, request):
    units = webutils.post(request, 'units', '')
    
    numcols = 0
    for field in request.POST:
        m = re.match('column_([0-9]+)_(type|label)', field)
        if m:
            col = int(m.group(1))
            if col + 1 > numcols:
                numcols = col + 1
                
    errors = []
    columns = {}
    
    data_labels = set()
    stddev_labels = set()
    stddev_cols = []
    for c in xrange(0, numcols):
        col = {'type':'', 'label':''}
        columns[c] = col
        col['type'] = webutils.post(request,'column_%d_type' % (c), "").strip()
        
        if col['type'] == "":
            errors.append(strings.experiment_upload_error_column_type_not_defined % (c+1))
            
        col['label'] = webutils.post(request, 'column_%d_label' % (c), "").strip()

        if col['label'] == "" and col['type'] in set(['stddev','data']):        
            errors.append(strings.experiment_upload_error_data_column_empty_label % (c+1))
            
        if col['label'] != "":
            if col['type'] == 'data':
                if col['label'] in data_labels:
                    errors.append(strings.experiment_upload_error_data_column_label_duplicated % (c+1))
                data_labels.add(col['label'])
            elif col['type'] == 'stddev':
                if col['label'] in stddev_labels:
                    errors.append(strings.experiment_upload_error_data_column_label_duplicated % (c+1))
                stddev_labels.add(col['label'])
                stddev_cols.append(c)
            else:
                col['label'] = ""

    [ errors.append(strings.experiment_upload_error_standard_deviation_label_does_not_match_any_data_column % (c+1, columns[c]['label'])) 
            for c in stddev_cols if columns[c]['label'] not in data_labels ]

    session.columns = []
    for c in xrange(0, numcols):
        col = upload.SessionColumn()
        col.type = columns[c]['type']
        col.label = columns[c]['label']
        col.column_number = c
        session.columns.append(col)
        
    session.units = units

    return {'columns':columns,'units':units}, [ uploadutils.ColumnError(e) for e in errors ]


@view_config(route_name='upload_config', renderer='ptmscout:/templates/upload/upload_config.pt')
def upload_config(request):
    session_id = int(request.matchdict['id'])
    session = upload.getSessionById(session_id, request.user) 
    submitted = webutils.post(request, 'submitted', "false") == "true"
    force = webutils.post(request, 'override', "false") != "false"
    
    allowoverride = False
    
    column_defs = []
    errors = []
    if submitted:
        commit = False
        column_defs, errors = parse_user_input(session, request)
        
        try:
            if len(errors) > 0:
                raise uploadutils.ErrorList(errors, True)
                
            uploadutils.check_data_column_assignments(session)
            commit = True
        except uploadutils.ErrorList, ce:
            allowoverride = not ce.critical
            commit = force and allowoverride
            errors = ce.error_list()
            
        if commit:
            session.stage = 'metadata'
            session.save()
            return HTTPFound(request.application_url + "/upload/%d/metadata" % (session_id))
    else:
        column_defs = uploadutils.assign_column_defaults(session)
    
    headers, data_rows = uploadutils.load_header_and_data_rows(session.data_file, uploadutils.MAX_ROW_CHECK)
    
    return {'allowoverride': allowoverride,
            'headers': headers,
            'data_rows': data_rows,
            'error': errors,
            'data_definitions':column_defs,
            'session_id': session_id,
            'column_values':upload.SessionColumn.column_values,
            'pageTitle': strings.experiment_upload_configure_page_title,
            'instruction': strings.experiment_upload_configure_message}
