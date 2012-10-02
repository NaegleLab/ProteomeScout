from mock import patch
from ptmscout.config import strings
from ptmscout.views.protein.data_view import protein_experiment_data_view
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockProtein, createMockUser, \
    createMockMeasurement, createMockExperiment, createMockPhosphopep, \
    createMockData
import unittest
from pyramid import testing

class TestProteinDataViews(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
        
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByProtein')
    @patch('ptmscout.database.protein.getProteinById')
    def test_protein_data_view_should_not_include_experiments_without_data(self, patch_getProtein, patch_getMods):
        request = DummyRequest()

        mock_prot = createMockProtein()
        request.matchdict['id'] = str(mock_prot.id)
        patch_getProtein.return_value = mock_prot
        
        mock_user = createMockUser("username", "email", "password", 1)
        request.user = mock_user
        
        mock_mod = createMockMeasurement(mock_prot.id, 24)
        exp = createMockExperiment(2, 0, 0)
        mock_mod.experiment = exp
        
        mock_mod2 = createMockMeasurement(mock_prot.id, 25)
        exp2 = createMockExperiment(3, 0, 0)
        mock_mod2.experiment = exp2
        
        pep1 = createMockPhosphopep(mock_prot.id)
        pep2 = createMockPhosphopep(mock_prot.id)
        pep3 = createMockPhosphopep(mock_prot.id)
        mock_mod.phosphopeps.append(pep1)
        mock_mod2.phosphopeps.append(pep2)
        mock_mod2.phosphopeps.append(pep3)

        mod_list = [mock_mod, mock_mod2]
        patch_getMods.return_value = mod_list 
        
        result = protein_experiment_data_view(request)
        
        patch_getProtein.assert_called_with(mock_prot.id)
        
        self.assertEqual(strings.protein_data_page_title, result['pageTitle'])
        self.assertEqual(mock_prot, result['protein'])
        
        self.assertEqual([], result['experiment_data'])
    
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByProtein')
    @patch('ptmscout.database.protein.getProteinById')
    def test_protein_data_view_get_experiment_data_for_protein(self, patch_getProtein, patch_getMods):
        request = DummyRequest()

        mock_prot = createMockProtein()
        request.matchdict['id'] = str(mock_prot.id)
        patch_getProtein.return_value = mock_prot
        
        mock_user = createMockUser("username", "email", "password", 1)
        request.user = mock_user
        
        mock_mod = createMockMeasurement(mock_prot.id, 24)
        exp = createMockExperiment(2, 0, 0)
        mock_mod.experiment = exp
        
        mock_mod2 = createMockMeasurement(mock_prot.id, 25)
        exp2 = createMockExperiment(3, 0, 0)
        mock_mod2.experiment = exp2
        
        pep1 = createMockPhosphopep(mock_prot.id)
        pep2 = createMockPhosphopep(mock_prot.id)
        pep3 = createMockPhosphopep(mock_prot.id)
        mock_mod.phosphopeps.append(pep1)
        mock_mod2.phosphopeps.append(pep2)
        mock_mod2.phosphopeps.append(pep3)
        
        data1 = createMockData(5, 'run1_12H', mock_mod.id)
        data2 = createMockData(5, 'run2_24H', mock_mod.id)
        mock_mod.data.extend(data1)
        mock_mod.data.extend(data2)

        data3 = createMockData(4, '24H_stuff', mock_mod2.id)
        mock_mod2.data.extend(data3)

        mod_list = [mock_mod, mock_mod2]
        patch_getMods.return_value = mod_list 
        
        result = protein_experiment_data_view(request)
        
        patch_getProtein.assert_called_with(mock_prot.id)
        
        self.assertEqual(strings.protein_data_page_title, result['pageTitle'])
        self.assertEqual(mock_prot, result['protein'])
        
        data1_s_vals = [(row.label, str(row.value), row.type) for row in sorted(data1, key=lambda item: item.priority)]        
        data2_s_vals = [(row.label, str(row.value), row.type) for row in sorted(data2, key=lambda item: item.priority)]
        data3_s_vals = [(row.label, str(row.value), row.type) for row in sorted(data3, key=lambda item: item.priority)]
        
        exp_data1 = [{'run':"run1_12H", 'phosphopeps': [pep1.getName()], 'values': data1_s_vals},
                     {'run':"run2_24H", 'phosphopeps': [pep1.getName()], 'values': data2_s_vals}]
        exp_exp1 = {'id': exp.id, 'title': exp.name, 'data':exp_data1}
        
        exp_data2 = [{'run':"24H_stuff", 'phosphopeps': [pep2.getName(), pep3.getName()], 'values': data3_s_vals}]
        exp_exp2 = {'id': exp2.id, 'title': exp2.name, 'data':exp_data2}
        
        self.assertEqual([exp_exp1, exp_exp2], result['experiment_data'])