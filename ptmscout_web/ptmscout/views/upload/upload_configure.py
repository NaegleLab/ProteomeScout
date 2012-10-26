from ptmscout.config import strings
from pyramid.view import view_config
from ptmscout.database import upload


def parse_user_input(session, request):
    pass


def get_columns_of_type(session, tp):
    cols = []
    for col in session.columns:
        if tp == col.type:
            cols.append(col)
    return cols

def check_data_rows(session):
    errors = []
    
    
    
    errors

class ColumnError(Exception):
    def __init__(self, message):
        self.message = message
        
    def __repr__(self):
        return self.message

def check_unique_column(session, ctype, required=False):
    cols = get_columns_of_type(session, ctype)

    if required and len(cols) == 0:
        raise ColumnError(strings.experiment_upload_warning_no_column_assignment % ctype)
    if len(cols) > 1:
        raise ColumnError(strings.experiment_upload_error_limit_one_column_of_type % ctype)
    
    if len(cols) == 0:
        return None
    return cols[0]

def check_stddev_maps_to_data(data_cols, stddev_cols):
    for col in stddev_cols:
        maps = False
        for c2 in data_cols:
            maps |= col.label == c2.label
        if not maps:
            raise ColumnError(strings.experiment_upload_error_standard_deviation_label_does_not_match_any_data_column % (col.label))
        
    
def check_data_column_assignments(session):
    acc_col = check_unique_column(session, 'accession', required=True)
    pep_col = check_unique_column(session, 'peptide', required=True)
    mod_col = check_unique_column(session, 'modification')
    species_col = check_unique_column(session, 'species')
    run_col = check_unique_column(session, 'run')
    
    data_cols   = get_columns_of_type(session, 'data')
    stddev_cols = get_columns_of_type(session, 'stddev')
    
    check_stddev_maps_to_data(data_cols, stddev_cols)
    
    
    

def assign_columns_by_name(session, header):
    pass

def assign_columns_from_session(session):
    pass

def assign_columns_from_session_history(header, session, current_user):
    pass

def assign_column_defaults(session, current_user):
    pass





def load_header_and_data_rows(session, N=20):
    pass




@view_config(route_name='upload_config', renderer='ptmscout:/templates/upload/upload_config.pt')
def upload_config(request):
    session_id = int(request.matchdict['id'])
    session = upload.getSessionById(session_id, request.user) 
    
    return {'pageTitle': strings.experiment_upload_configure_page_title,
            'instruction': strings.experiment_upload_configure_message}