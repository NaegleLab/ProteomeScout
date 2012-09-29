from pyramid import testing
import unittest
from pyramid.testing import DummyRequest
from ptmscout.experiments import experiment_listing, view_experiment,\
    browse_experiment
from mock import patch, Mock
from ptmscout import strings
from tests.views.test_user_management import createUserForTest
from pyramid.httpexceptions import HTTPForbidden
from tests.views.mocking import createMockExperiment, createMockPermission,\
    createMockProtein, createMockModification, createMockPhosphopep
from ptmscout.database.protein import Species

class ExperimentViewTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()
    
    @patch('ptmscout.database.experiment.getExperimentTree')
    def test_experiment_listing(self, patch_getExperimentTree):
        request = DummyRequest()
        exp1 = createMockExperiment(1,1)
        exp2 = createMockExperiment(2,1)
        
        patch_getExperimentTree.return_value = [exp1,exp2]
        request.user = None
        
        result = experiment_listing(request)
        
        patch_getExperimentTree.assert_called_once_with(request.user)
        self.assertEqual([exp1,exp2], result['experiments'])
        self.assertEqual(strings.experiments_page_title, result['pageTitle'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_view(self, patch_getExperiment):
        mock_experiment = Mock('ptmscout.database.experiment.Experiment')
        patch_getExperiment.return_value = mock_experiment
        
        mock_experiment.name = "experiment name"
        
        request = DummyRequest()
        request.user = None
        request.matchdict['id'] = 1 
        
        parameters = view_experiment(request)
        
        patch_getExperiment.assert_called_once_with(1, request.user)
        self.assertEqual(mock_experiment, parameters['experiment'])
        self.assertEqual(strings.experiment_page_title % ("experiment name"), parameters['pageTitle'])
    
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_browse_view_should_show_error_when_search_acc_gene_empty(self, patch_getExperiment):
        mock_experiment = createMockExperiment(1, 0, 0)
        patch_getExperiment.return_value = mock_experiment
        
        mock_experiment.name = "experiment name"
        
        request = DummyRequest()
        
        request.POST['submitted'] = "true"
        request.POST['acc_search'] = "   "
        request.POST['stringency'] = "2 "

        request.user = None
        request.matchdict['id'] = 1
    
        result = browse_experiment(request)
    
        self.assertEqual([], result['proteins'])
        self.assertEqual({}, result['modifications'])
        self.assertEqual(True, result['submitted'])
        
        self.assertEqual(mock_experiment, result['experiment'])
        self.assertEqual(strings.experiment_browse_page_title % ("experiment name"), result['pageTitle'])
        
        self.assertEqual("", result['acc_search'])
        self.assertEqual("2", result['stringency'])
        self.assertEqual(mock_experiment, result['experiment'])
    
    @patch('ptmscout.database.modifications.getModificationsByExperiment')
    @patch('ptmscout.database.protein.getProteinsByAccession')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_browse_view_when_searching(self, patch_getExperiment, patch_getProteins, patch_getMods):
        mock_experiment = createMockExperiment(1, 0, 0)
        patch_getExperiment.return_value = mock_experiment
        
        print "test search..."
        mock_experiment.name = "experiment name"
        
        request = DummyRequest()
        request.POST['submitted'] = "true"
        request.POST['acc_search'] = "ACK1"
        request.POST['stringency'] = "2"
        request.matchdict['id'] = 1
        
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        
        p1 = createMockProtein()
        protein_list = [p1]
        
        patch_getProteins.return_value = protein_list
        
        m1 = createMockModification(p1.id, 1)
        m2 = createMockModification(p1.id, 4)
        
        m1.protein = p1
        m2.protein = p1
        
        m1.experiment_id = mock_experiment.id
        m2.experiment_id = mock_experiment.id
        
        pep1 = createMockPhosphopep(p1.id)
        pep2 = createMockPhosphopep(p1.id)
        pep3 = createMockPhosphopep(p1.id)
        
        m1.phosphopeps.append(pep1)
        m1.phosphopeps.append(pep2)
        m2.phosphopeps.append(pep3)
        
        mod_list = [m1, m2]
        
        patch_getMods.return_value = mod_list
        
        result = browse_experiment(request)
        
        patch_getProteins.assert_called_once_with(["ACK1"])
        patch_getExperiment.assert_called_once_with(1, request.user)
        patch_getMods.assert_called_once_with(1, request.user, [p1.id])
        self.assertEqual("ACK1", result['acc_search'])
        self.assertEqual("2", result['stringency'])
        self.assertEqual(mock_experiment, result['experiment'])
        self.assertEqual(strings.experiment_browse_page_title % ("experiment name"), result['pageTitle'])
        
        sorted_peps = sorted([pep1, pep2, pep3], key=lambda pep: pep.getName())
        parsed_peps = [{'site':pep.getName(), 'peptide':pep.getPeptide()} for pep in sorted_peps]
        
        self.assertEqual(protein_list, result['proteins'])
        self.assertEqual({p1.id: parsed_peps}, result['modifications'])
        self.assertEqual(True, result['submitted'])
        
    
    @patch('ptmscout.database.modifications.getModificationsByExperiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_experiment_browse_view_when_no_search_parameters(self, patch_getExperiment, patch_getMods):
        mock_experiment = createMockExperiment(1, 0, 0)
        patch_getExperiment.return_value = mock_experiment
        
        mock_experiment.name = "experiment name"
        
        request = DummyRequest()
        p1 = createMockProtein()
        ptm_user = createUserForTest("username", "email", "password", 1)
        request.user = ptm_user
        request.matchdict['id'] = 1
        protein_list = [p1]
        
        m1 = createMockModification(p1.id, 1)
        m2 = createMockModification(p1.id, 4)
        
        m1.protein = p1
        m2.protein = p1
        
        m1.experiment_id = mock_experiment.id
        m2.experiment_id = mock_experiment.id
        
        pep1 = createMockPhosphopep(p1.id)
        pep2 = createMockPhosphopep(p1.id)
        pep3 = createMockPhosphopep(p1.id)
        
        m1.phosphopeps.append(pep1)
        m1.phosphopeps.append(pep2)
        m2.phosphopeps.append(pep3)
        
        mod_list = [m1, m2]
        

        
        patch_getMods.return_value = mod_list
        
        result = browse_experiment(request)
        
        patch_getExperiment.assert_called_once_with(1, request.user)
        patch_getMods.assert_called_once_with(1, request.user)
        self.assertEqual("", result['acc_search'])
        self.assertEqual("1", result['stringency'])
        self.assertEqual(mock_experiment, result['experiment'])
        self.assertEqual(strings.experiment_browse_page_title % ("experiment name"), result['pageTitle'])
        
        sorted_peps = sorted([pep1, pep2, pep3], key=lambda pep: pep.getName())
        parsed_peps = [{'site':pep.getName(), 'peptide':pep.getPeptide()} for pep in sorted_peps]
        
        self.assertEqual(protein_list, result['proteins'])
        self.assertEqual({p1.id: parsed_peps}, result['modifications'])
        self.assertEqual(False, result['submitted'])