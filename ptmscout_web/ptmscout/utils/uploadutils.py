from ptmscout.config import strings, settings
import csv
import os
from ptmscout.utils.webutils import call_catch
import re
import sys
from ptmscout.utils import protein_utils

MAX_ROW_CHECK=100


class ErrorList(Exception):
    def __init__(self, errors, critical=True):
        self.errors = errors
        self.critical = critical
        
    def __repr__(self):
        st = ""
        for err in self.errors:
            st += "%s\n" % (err)
    
    def error_list(self):
        return [ e.message for e in self.errors ]
            
class ParseError(Exception):
    def __init__(self, row, col, msg):
        self.row = row
        self.col = col
        self.msg = msg
    
    def __repr__(self):
        if(self.col != None):
            return "Line %d, Column %d: %s" % (self.row, self.col, self.msg)
        else:
            return "Line %d: %s" % (self.row, self.msg)
    
    message = property(__repr__)    

class ColumnError(Exception):
    def __init__(self, message):
        self.message = message
        
    def __repr__(self):
        return self.message

def get_columns_of_type(session, tp):
    cols = []
    for col in session.columns:
        if tp == col.type:
            cols.append(col)
    return cols



def check_unique_column(session, ctype, required=False):
    cols = get_columns_of_type(session, ctype)

    if required and len(cols) == 0:
        raise ColumnError(strings.experiment_upload_error_no_column_assignment % (ctype,))
    if len(cols) > 1:
        raise ColumnError(strings.experiment_upload_error_limit_one_column_of_type % (ctype,))
    
    if len(cols) == 0:
        return None
    return cols[0]


def check_data_rows(session, acc_col, pep_col, mod_col, run_col, data_cols, stddev_cols, N=MAX_ROW_CHECK):
    errors = []
    header, data = load_header_and_data_rows(session, N)
    
    keys = set([])
    
    r = 0
    for row in data:
        if len(row) < len(header):
            continue
        
        r+=1
        accession = row[acc_col.column_number].strip()
        peptide = row[pep_col.column_number].strip()
        modification = row[mod_col.column_number].strip()

        acc_type = protein_utils.get_accession_type(accession)
        if acc_type not in protein_utils.get_valid_accession_types():
            errors.append(ParseError(r, acc_col.column_number+1, strings.experiment_upload_warning_acc_column_contains_bad_accessions))
                
        if not protein_utils.check_peptide_alphabet(peptide):
            errors.append(ParseError(r, pep_col.column_number+1, strings.experiment_upload_warning_peptide_column_contains_bad_peptide_strings))
            
        if not protein_utils.check_peptide_modification_valid(peptide, modification):
            errors.append(ParseError(r, mod_col.column_number+1, strings.experiment_upload_warning_modifications_do_not_match_amino_acids))
        
        run = None
        if run_col != None:
            run = row[run_col.column_number].strip()
            k = (accession, peptide, run)
            if k in keys:
                errors.append(ParseError(r, None, strings.experiment_upload_warning_full_dupe))
            keys.add(k)
        else:
            k = (accession, peptide)
            if k in keys:
                errors.append(ParseError(r, None, strings.experiment_upload_warning_no_run_column))
            keys.add(k)
            
        for c in data_cols + stddev_cols:
            try:
                float(row[c.column_number].strip())
            except:
                errors.append(ParseError(r, c.column_number+1, strings.experiment_upload_warning_data_column_not_numeric))
        
        
    
    return errors


def check_data_column_assignments(session):
    errors = []
    
    acc_col     = call_catch(ColumnError, errors, check_unique_column, session, 'accession', required=True)
    pep_col     = call_catch(ColumnError, errors, check_unique_column, session, 'peptide', required=True)
    mod_col     = call_catch(ColumnError, errors, check_unique_column, session, 'modification', required=True)
    run_col     = call_catch(ColumnError, errors, check_unique_column, session, 'run')
    
    critical = True
    
    if len(errors) == 0:
        critical = False
        data_cols   = get_columns_of_type(session, 'data')
        stddev_cols = get_columns_of_type(session, 'stddev')
        errors.extend(check_data_rows(session, acc_col, pep_col, mod_col, run_col, data_cols, stddev_cols))
        
    if len(errors) > 0:
        raise ErrorList(errors, critical)
    


def get_column_from_header(header):
    for h in header:
        h = h.lower()
        col = {'type':'','label':''}
        
        m = re.match('^(data|stddev):(.+):(.+)$', h)
        if m:
            col['type'] = 'stddev' if m.group(1) == 'stddev' or m.group(2).find('stddev') == 0 else 'data' 
            col['label'] = m.group(3)
        elif h.find('data') == 0:
            col['type'] = 'data'
        elif h.find('stddev') == 0:
            col['type'] = 'stddev'
        elif h.find('acc') == 0:
            col['type'] = 'accession'
        elif h.find('mod') == 0:
            col['type'] = 'modification'
        elif h.find('pep') == 0:
            col['type'] = 'peptide'
        elif h == 'run':
            col['type'] = 'run'
        else:
            col['type'] = 'none'
        
        
        yield col
    
def find_units(header):
    for h in header:
        m = re.match('^data:(.+):.+$',h)
        if m:
            return m.group(1)
    return ""

def assign_columns_by_name(header):
    columns = {}

    c = 0
    for col in get_column_from_header(header):
        columns[c] = col
        c += 1
    
    units = find_units(header)
    
    return {'columns':columns, 'units':units}

def assign_columns_from_session(session):
    columns = {}
    
    for col in session.columns:
        columns[col.column_number] = {'type':col.type, 'label':col.label}

    return {'columns':columns, 'units':session.units}

def assign_columns_from_session_history(session, header):
    history_session = session.getAncestor()
    result = assign_columns_from_session(history_session)
    
    keys = result['columns'].keys()
    for c in keys:
        if not c < len(header):
            del result['columns'][c]
    
    return result

def assign_column_defaults(session):
    header, _ = load_header_and_data_rows(session)
    
    if session.columns != []:
        return assign_columns_from_session(session)
    
    if session.parent_experiment != None:
        return assign_columns_from_session_history(session, header)
    
    return assign_columns_by_name(header)


def load_header_and_data_rows(session, N=-1):
    ifile = csv.reader(open(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, session.data_file), 'rb'), delimiter='\t')
    i = 0
    
    header = ifile.next()
    
    width = len(header)
    while(width > 0 and header[width-1].strip() == ''):
        width-=1

    header = header[0:width]
    
    rows = []
    for row in ifile:
        if i >= N:
            break
        row = row[0:width]
        rows.append(row)
    
    return header, rows
