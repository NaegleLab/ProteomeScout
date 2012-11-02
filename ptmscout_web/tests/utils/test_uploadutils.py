import unittest
from ptmscout.utils.uploadutils import ColumnError,\
    check_data_column_assignments, ErrorList, assign_column_defaults,\
    assign_columns_by_name, assign_columns_from_session_history, check_data_rows,\
    ParseError, check_modification_type_matches_peptide \
    
from tests.views.mocking import createMockUser, createMockSession,\
    createMockSessionColumn, createMockPTM
from mock import patch, call
from ptmscout.config import strings

class TestUploadUtils(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)

    @patch('ptmscout.database.modifications.findMatchingPTM')
    def test_check_modification_type_matches_peptide_should_throw_error_when_mismatch_of_number_of_mods(self, patch_findPTM):
        mods = []
        patch_findPTM.return_value = mods, True
        
        try:
            check_modification_type_matches_peptide(1, "DFERtGdYDFqERAS", "Sulfation, METH")
        except ParseError, pe:
            self.assertEqual(1, pe.row)
            self.assertEqual(strings.experiment_upload_warning_wrong_number_of_mods % (2, 3), pe.msg)
        else:
            self.fail("Expected exception ParseError")
        
    @patch('ptmscout.database.modifications.findMatchingPTM')
    def test_check_modification_type_matches_peptide_should_throw_error_when_no_mod_of_type(self, patch_findPTM):
        mods = []
        patch_findPTM.return_value = mods, False
        
        try:
            check_modification_type_matches_peptide(1, "DFERTGDYDFqERAS", "Sulfation")
        except ParseError, pe:
            self.assertEqual(1, pe.row)
            self.assertEqual(strings.experiment_upload_warning_modifications_not_valid % ('Sulfation'), pe.msg)
        else:
            self.fail("Expected exception ParseError")
            
        patch_findPTM.assert_any_call("Sulfation", "Q") 

    @patch('ptmscout.database.modifications.findMatchingPTM')
    def test_check_modification_type_matches_peptide_should_throw_error_when_no_matching_mod(self, patch_findPTM):
        mods = []
        patch_findPTM.return_value = mods, True
        
        try:
            check_modification_type_matches_peptide(1, "DFERTGDYDFqERAS", "Sulfation")
        except ParseError, pe:
            self.assertEqual(1, pe.row)
            self.assertEqual(strings.experiment_upload_warning_modifications_do_not_match_amino_acids % ('Sulfation', 'Q'), pe.msg)
        else:
            self.fail("Expected exception ParseError")
            
        patch_findPTM.assert_any_call("Sulfation", "Q") 

    @patch('ptmscout.database.modifications.findMatchingPTM')
    def test_check_modification_type_matches_peptide_should_not_throw_error_on_ambiguous_modification_if_category_available(self, patch_findPTM):
        ptm1 = createMockPTM(name="Methylation", keywords=["METH"])
        ptm2 = createMockPTM(name="2-methylglutamine", target="Q", parent=ptm1, keywords=["methylation", "METH"])
        ptm3 = createMockPTM(name="N5-methylglutamine", target="Q", parent=ptm1, keywords=["methylation", "METH"])
        
        mods = [ptm1,ptm2,ptm3]
        patch_findPTM.return_value = mods, True
        
        check_modification_type_matches_peptide(1, "DFERTGDYDFqERAS", "METH")
        
        patch_findPTM.assert_any_call("METH", "Q") 
    
    @patch('ptmscout.database.modifications.findMatchingPTM')
    def test_check_modification_type_matches_peptide_should_throw_error_on_ambiguous_modification_if_no_fallback_available(self, patch_findPTM):
        ptm1 = createMockPTM(name="2-methylglutamine", target="Q", keywords=["methylation", "METH"])
        ptm2 = createMockPTM(name="N5-methylglutamine", target="Q", keywords=["methylation", "METH"])
        
        mods = [ptm1,ptm2]
        patch_findPTM.return_value = mods, True
        
        try:
            check_modification_type_matches_peptide(1, "DFERTGDYDFqERAS", "METH")
        except ParseError, pe:
            self.assertEqual(1, pe.row)
            self.assertEqual(strings.experiment_upload_warning_ambiguous_modification_type_for_amino_acid % ("METH", "Q"), pe.msg) 
        else:
            self.fail("Expected exception ParseError")
        
        patch_findPTM.assert_any_call("METH", "Q") 
    
    @patch('ptmscout.database.modifications.findMatchingPTM')
    def test_check_modification_type_matches_should_succeed(self, patch_findPTM):
        ptm1 = createMockPTM(name="Phosphorylation", keywords=["PHOS"])
        ptm2 = createMockPTM(name="phosphoserine", target="S", parent=ptm1, keywords=["Phosphorylation", "PHOS"])
        
        patch_findPTM.return_value = [ptm1,ptm2], True
               
        check_modification_type_matches_peptide(1, "AKsPVPKsPVEEK", "PHOS")
        
        patch_findPTM.assert_any_call("PHOS", "S")

    @patch('ptmscout.utils.uploadutils.check_modification_type_matches_peptide')
    def test_check_data_rows_should_compile_list_of_errors(self, patch_test_peptide):
        user = createMockUser()
        session = createMockSession(user)
        
        patch_test_peptide.side_effect = ParseError(1, None, "Peptide modification error")
        
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
        
        expected_errors = ["Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 2, Column 8: " + strings.experiment_upload_warning_data_column_not_numeric,
                           "Line 1: Peptide modification error",
                           "Line 4, Column 1: " + strings.experiment_upload_warning_acc_column_contains_bad_accessions,
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 9, Column 4: " + strings.experiment_upload_warning_peptide_column_contains_bad_peptide_strings,
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 12: " + strings.experiment_upload_warning_no_run_column,
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error",
                           "Line 1: Peptide modification error"]
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
        headers = ['Acc', 'pep', 'other', 'gene(name)', \
                   'run', 'moD_type', 'modification', \
                   'data', 'data:time(min):20', \
                   'data:stddev(sec):10', 'stddev', 'stddev:time(min):20']
        
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
                                     11:{'type':'stddev','label':'20'},
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
    def test_check_data_column_assignments_should_raise_errors(self, patch_check_unique, patch_check_rows):
        ce = ColumnError("Some error")
        patch_check_unique.return_value=createMockSessionColumn(0, 'accession', 1)
        patch_check_rows.return_value = [ce]
        
        user = createMockUser()
        session = createMockSession(user)
        session.columns = []
        try:
            check_data_column_assignments(session)
        except ErrorList, el:
            self.assertEqual([ce.message], el.error_list())
            self.assertEqual(False, el.critical)
        else:
            self.fail("Expected exception: ErrorList([ColumnError()])")
            
    @patch('ptmscout.utils.uploadutils.check_data_rows')
    @patch('ptmscout.utils.uploadutils.check_unique_column')
    def test_check_data_column_assignments_should_raise_error_not_check_rows(self, patch_check_unique, patch_check):
        ce = ColumnError("Some error")
        patch_check_unique.side_effect = ce
        patch_check.return_value = []
        
        user = createMockUser()
        session = createMockSession(user)
        session.columns = []
        
        try:
            check_data_column_assignments(session)
        except ErrorList, el:
            self.assertEqual([ce.message] * 4, el.error_list())
            self.assertEqual(True, el.critical)
        else:
            self.fail("Expected exception: ErrorList([ColumnError()])")

        self.assertFalse(patch_check.called)