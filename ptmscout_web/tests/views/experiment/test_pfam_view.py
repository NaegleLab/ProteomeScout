from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
from mock import patch
from ptmscout.views.experiment.pfam_view import show_pfam_view,\
    format_pfam_domains, format_pfam_sites
from tests.views.mocking import createMockExperiment, createMockProtein,\
    createMockMeasurement, createMockPeptide, createMockDomain,\
    createMockPeptideModification, createMockPTM, createMockUser
from ptmscout.config import strings
import base64
import json

class TestPfamView(UnitTestCase):
    
    def create_test_data(self):
        p1 = createMockProtein()
        p2 = createMockProtein()
        p3 = createMockProtein()
        
        d1 = createMockDomain(p1.id)
        d2 = createMockDomain(p1.id)
        d3 = createMockDomain(p2.id)
        
        d1.label = 'pfam_domain_1'
        d2.label = 'pfam_domain_2'
        d3.label = 'pfam_domain_1'
        
        p1.domains = [d1,d2]
        p2.domains = [d3]
        p3.domains = []
        
        m1 = createMockMeasurement(p1.id, 1)
        m2 = createMockMeasurement(p1.id, 1)
        m3 = createMockMeasurement(p2.id, 1)
        m4 = createMockMeasurement(p3.id, 1)
        
        m1.protein = p1
        m2.protein = p1
        m3.protein = p2
        m4.protein = p3
        
        pep1 = createMockPeptide(p1.id)
        pep2 = createMockPeptide(p1.id)
        pep3 = createMockPeptide(p1.id)
        pep4 = createMockPeptide(p2.id)
        pep5 = createMockPeptide(p2.id)
        
        d1 = createMockDomain(p1.id, label='pfam001')
        d2 = createMockDomain(p1.id, label='pfam002')
        d3 = createMockDomain(p1.id, label='pfam004')
        
        pep1.protein_domain = d1
        pep2.protein_domain = d2
        pep3.protein_domain = d1
        pep4.protein_domain = d3
        pep5.protein_domain = None  
        
        mod = createMockPTM()
        
        createMockPeptideModification(m1, pep1, mod)
        createMockPeptideModification(m1, pep2, mod)
        createMockPeptideModification(m2, pep2, mod)
        createMockPeptideModification(m2, pep3, mod)
        createMockPeptideModification(m3, pep4, mod)
        createMockPeptideModification(m3, pep5, mod)
        
        measurements = [m1,m2,m3,m4]
        
        return measurements
    
    
    def test_format_pfam_domains(self):
        measurements = self.create_test_data()
        result = format_pfam_domains(measurements)
        
        expected_table = [('pfam_domain_1', 2), ('pfam_domain_2', 1), ('None', 1)] 
        json_table = base64.b64encode(json.dumps(expected_table))

        self.assertEqual({'table':expected_table, 'json':json_table}, result)
        
    def test_format_pfam_sites(self):
        measurements = self.create_test_data()
        result = format_pfam_sites(measurements)
        
        expected_table = [('pfam002', 2), ('pfam001', 2), ('None', 1), ('pfam004', 1)]
        json_table = base64.b64encode(json.dumps(expected_table))
        
        self.assertEqual({'table':expected_table, 'json':json_table}, result)
            
    
    
    @patch('ptmscout.views.experiment.pfam_view.format_pfam_domains')
    @patch('ptmscout.views.experiment.pfam_view.format_pfam_sites')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_view(self, patch_getExperiment, patch_getMeasurements, patch_formatSites, patch_formatDomains):
        request = DummyRequest()
        expid = 1
        exp = createMockExperiment(expid, 0, 0)
        request.matchdict['id'] = expid
        request.user = createMockUser()
        
        patch_getExperiment.return_value = exp
        
        measurements = ["some", "measurements"]
        patch_getMeasurements.return_value = measurements
        
        sites = ["some", "formatted", "sites"]
        patch_formatSites.return_value = sites
        domains = ["some","formatted","domains"]
        patch_formatDomains.return_value = domains
        
        request.user.experimentOwner.return_value = True
        
        result = show_pfam_view(request)
        
        request.user.experimentOwner.assert_called_once_with(exp)
        
        patch_getExperiment.assert_called_once_with(expid, request.user)
        patch_getMeasurements.assert_called_once_with(expid, request.user)
        
        patch_formatSites.assert_called_once_with(measurements)
        patch_formatDomains.assert_called_once_with(measurements)
        
        self.assertEqual(strings.experiment_pfam_page_title % (exp.name), result['pageTitle'])
        self.assertEqual(exp, result['experiment'])
        
        self.assertEqual(True, result['user_owner'])
        self.assertEqual(sites, result['sites'])
        self.assertEqual(domains, result['domains'])
        
        
class IntegrationTestPfamView(IntegrationTestCase):
    def test_view_integration(self):
        self.ptmscoutapp.get("/experiments/26/pfam", status=200)