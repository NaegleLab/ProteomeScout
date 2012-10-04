from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from pyramid.testing import DummyRequest
from mock import patch
from ptmscout.views.experiment.pfam_view import show_pfam_view,\
    format_pfam_domains, format_pfam_sites
from tests.views.mocking import createMockExperiment, createMockProtein,\
    createMockMeasurement, createMockPhosphopep, createMockDomain
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
        d4 = createMockDomain(p2.id)
        
        d1.label = 'pfam_domain_1'
        d2.label = 'pfam_domain_2'
        d3.label = 'pfam_domain_1'
        d4.label = '~~~'
        
        p1.domains = [d1,d2,d4]
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
        
        pep1 = createMockPhosphopep(p1.id)
        pep2 = createMockPhosphopep(p1.id)
        pep3 = createMockPhosphopep(p1.id)
        pep4 = createMockPhosphopep(p2.id)
        pep5 = createMockPhosphopep(p2.id)
        
        pep1.pfam_site = 'pfam001'
        pep2.pfam_site = 'pfam002'
        pep3.pfam_site = 'pfam001'
        pep4.pfam_site = 'pfam004'
        pep5.pfam_site = '~~~'  
        
        m1.phosphopeps = [pep1, pep2]
        m2.phosphopeps = [pep2, pep3]
        m3.phosphopeps = [pep4, pep5]
        
        measurements = [m1,m2,m3]
        
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
        request.user = None
        
        patch_getExperiment.return_value = exp
        
        measurements = ["some", "measurements"]
        patch_getMeasurements.return_value = measurements
        
        sites = ["some", "formatted", "sites"]
        patch_formatSites.return_value = sites
        domains = ["some","formatted","domains"]
        patch_formatDomains.return_value = domains
        
        result = show_pfam_view(request)
        
        patch_getExperiment.assert_called_once_with(expid, request.user)
        patch_getMeasurements.assert_called_once_with(expid, request.user)
        
        patch_formatSites.assert_called_once_with(measurements)
        patch_formatDomains.assert_called_once_with(measurements)
        
        self.assertEqual(strings.experiment_pfam_page_title % (exp.name), result['pageTitle'])
        self.assertEqual(exp, result['experiment'])
        
        self.assertEqual(sites, result['sites'])
        self.assertEqual(domains, result['domains'])
        
        
class IntegrationTestPfamView(IntegrationTestCase):
    def test_view_integration(self):
        self.ptmscoutapp.get("/experiments/26/pfam", status=200)