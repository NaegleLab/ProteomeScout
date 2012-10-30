import unittest
from ptmscout.utils.uploadutils import ColumnError,\
    check_data_column_assignments, ErrorList, assign_column_defaults,\
    assign_columns_by_name, assign_columns_from_session_history, check_data_rows
from tests.views.mocking import createMockUser, createMockSession,\
    createMockSessionColumn
from mock import patch
from ptmscout.config import strings

class TestUploadUtils(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def test_check_data_rows_should_compile_list_of_errors(self):
        user = createMockUser()
        session = createMockSession(user)
        c1=createMockSessionColumn(0, 'accession', session.id)
        c2=createMockSessionColumn(1, 'none', session.id)
        c3=createMockSessionColumn(2, 'modification', session.id)
        c4=createMockSessionColumn(3, 'peptide', session.id)
        c5=createMockSessionColumn(4, 'data', session.id, label='0')
        c6=createMockSessionColumn(5, 'data', session.id, label='30')
        c7=createMockSessionColumn(6, 'data', session.id, label='60')
        c8=createMockSessionColumn(7, 'data', session.id, label='120')
        c9=createMockSessionColumn(8, 'data', session.id, label='240')
        c10=createMockSessionColumn(9, 'data', session.id, label='480')
        c11=createMockSessionColumn(10, 'data', session.id, label='960')
        c12=createMockSessionColumn(11, 'data', session.id, label='1440')
        c13=createMockSessionColumn(12, 'stddev', session.id, label='0')
        c14=createMockSessionColumn(13, 'stddev', session.id, label='30')
        c15=createMockSessionColumn(14, 'stddev', session.id, label='60')
        c16=createMockSessionColumn(15, 'stddev', session.id, label='120')
        c17=createMockSessionColumn(16, 'stddev', session.id, label='240')
        c18=createMockSessionColumn(17, 'stddev', session.id, label='480')
        c19=createMockSessionColumn(18, 'stddev', session.id, label='960')
        c20=createMockSessionColumn(19, 'stddev', session.id, label='1440')
        c21=createMockSessionColumn(20, 'none', session.id)
        c22=createMockSessionColumn(21, 'none', session.id)
        c23=createMockSessionColumn(22, 'none', session.id)
        
        session.columns = [c1,c2,c3,c4,c5,c6,c7,c8,c9,\
                           c10,c11,c12,c13,c14,c15,\
                           c16,c17,c18,c19,c20,c21,c22,c23]
        
        session.data_file = 'test/test_dataset_formatted_w_errors.txt'
        
        errors = check_data_rows(session, c1, c4, c3, None, [c5,c6,c7,c8,c9,c10,c11,c12], [c13,c14,c15,c16,c17,c18,c19,c20])
        
        expected_errors = [
                           "Line 2, Column 8: " + strings.experiment_upload_warning_data_column_not_numeric,
                           "Line 4, Column 1: " + strings.experiment_upload_warning_acc_column_contains_bad_accessions,
                           "Line 9, Column 4: " + strings.experiment_upload_warning_peptide_column_contains_bad_peptide_strings,
                           "Line 12: " + strings.experiment_upload_warning_no_run_column]
        self.maxDiff = None        
        
        self.assertEqual(expected_errors, [ e.message for e in errors ] )
    
    
    def test_assign_columns_from_session_history_should_truncate_if_too_many(self):
        headers = ['acc', 'pep', 'data', 'gene(name)']
        
        user = createMockUser()
        session = createMockSession(user)
        history_session = createMockSession(user)
        
        c1=createMockSessionColumn(0, 'accession', history_session.id)
        c2=createMockSessionColumn(1, 'peptide', history_session.id)
        c5=createMockSessionColumn(2, 'data', history_session.id, label='20')
        c3=createMockSessionColumn(3, 'none', history_session.id)
        c4=createMockSessionColumn(4, 'none', history_session.id)
        
        history_session.columns = [c1,c2,c3,c4,c5]
        history_session.units = 'historical_units'
        session.getAncestor.return_value = history_session
        
        defs = assign_columns_from_session_history(session, headers)
        
        session.getAncestor.assert_called_once_with()
        
        self.assertEqual({
                          'columns':{0:{'type':'accession','label':''},
                                     1:{'type':'peptide','label':''},
                                     2:{'type':'data','label':'20'},
                                     3:{'type':'none','label':''},
                                     },
                          'units':'historical_units'
                          }, defs)
        
        
        

    def test_assign_columns_by_name_should_generate_correct_columns(self):
        headers = ['acc', 'pep', 'other', 'gene(name)', \
                   'run', 'mod_type', 'modification', \
                   'data', 'data:time(min):20', \
                   'stddev:time(sec):10', 'stddev']
        
        defs = assign_columns_by_name(headers)
        
        self.assertEqual({
                          'columns':{0:{'type':'accession','label':''},
                                     1:{'type':'peptide','label':''},
                                     2:{'type':'none','label':''},
                                     3:{'type':'none','label':''},
                                     4:{'type':'run','label':''},
                                     5:{'type':'modification','label':''},
                                     6:{'type':'modification','label':''},
                                     7:{'type':'data','label':''},
                                     8:{'type':'data','label':'20'},
                                     9:{'type':'stddev','label':'10'},
                                     10:{'type':'stddev','label':''},
                                     },
                          'units':'time(min)'
                          }, defs)
        
        
        

    @patch('ptmscout.utils.uploadutils.load_header_and_data_rows')        
    @patch('ptmscout.utils.uploadutils.assign_columns_from_session_history')
    def test_assign_column_defaults_should_assign_from_history_if_session_empty(self, patch_assign_columns, patch_load_header):
        user = createMockUser()
        session = createMockSession(user)
        session.parent_experiment = 10
        session.columns = []
        
        header = "some header"
        patch_load_header.return_value = header, []
        
        patch_assign_columns.return_value = {"some assignments":"vals"}
        rval = assign_column_defaults(session)
        
        self.assertEqual({"some assignments":"vals"}, rval)
        patch_assign_columns.assert_called_with(session, header)
        patch_load_header.assert_called_once_with(session)


    @patch('ptmscout.utils.uploadutils.load_header_and_data_rows')        
    @patch('ptmscout.utils.uploadutils.assign_columns_from_session')
    def test_assign_column_defaults_should_assign_from_session_if_not_empty(self, patch_assign_columns, patch_load_header):
        user = createMockUser()
        session = createMockSession(user)
        session.parent_experiment = None
        session.columns = ["some columns"]
        
        header = "some header"
        patch_load_header.return_value = header, []
        
        patch_assign_columns.return_value = {"some assignments":"vals"}
        rval = assign_column_defaults(session)
        
        self.assertEqual({"some assignments":"vals"}, rval)
        patch_assign_columns.assert_called_with(session)
        patch_load_header.assert_called_once_with(session)


    @patch('ptmscout.utils.uploadutils.load_header_and_data_rows')        
    @patch('ptmscout.utils.uploadutils.assign_columns_by_name')
    def test_assign_column_defaults_should_assign_from_data_file_headers_if_session_empty_and_no_history(self, patch_assign_columns, patch_load_header):
        user = createMockUser()
        session = createMockSession(user)
        session.parent_experiment = None
        session.columns = []
        
        header = "some header"
        patch_load_header.return_value = header, []
        
        patch_assign_columns.return_value = {"some assignments":"vals"}
        rval = assign_column_defaults(session)
        
        self.assertEqual({"some assignments":"vals"}, rval)
        patch_assign_columns.assert_called_with(header)
        patch_load_header.assert_called_once_with(session)

    @patch('ptmscout.utils.uploadutils.check_data_rows')
    @patch('ptmscout.utils.uploadutils.check_unique_column')
    def test_check_data_column_assignments_should_raise_error(self, patch_check_unique, patch_check):
        ce = ColumnError("Some error")
        patch_check_unique.side_effect = ce
        patch_check.return_value = [ce]
        
        user = createMockUser()
        session = createMockSession(user)
        session.columns = []
        
        try:
            check_data_column_assignments(session)
        except ErrorList, el:
            self.assertEqual([ce.message] * 5, el.error_list())
        else:
            self.fail("Expected exception: ErrorList([ColumnError()])")

