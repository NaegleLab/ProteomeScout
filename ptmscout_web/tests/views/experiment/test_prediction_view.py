from pyramid import testing
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockProtein, createMockMeasurement,\
    createMockPeptide, createMockExperiment, createMockScansite, createMockPTM,\
    createMockPeptideModification, createMockUser
from ptmscout.views.experiment.prediction_view import prediction_view,\
    format_predictions, filter_predictions
from ptmscout.config import strings
from mock import patch
import base64
import json
from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase

class TestPredictionViews(UnitTestCase):
    def setUp(self):
        self.config = testing.setUp()
        
    def tearDown(self):
        testing.tearDown()
    
    def test_filter_predictions(self):
        pred1 = createMockScansite(1, 1)
        pred2 = createMockScansite(1, 1)
        
        pred1.percentile = 3
        pred2.percentile = 0.2
        
        result = filter_predictions([pred1,pred2])
        
        self.assertEqual([pred2], result)
        

    def createMockMeasurements(self, expid):
        exp = createMockExperiment(expid, 0, 0)
        
        p1 = createMockProtein()
        m1 = createMockMeasurement(p1.id, expid)
        m2 = createMockMeasurement(p1.id, expid)
        m3 = createMockMeasurement(p1.id, expid)
        
        pep1 = createMockPeptide(p1.id)
        pep2 = createMockPeptide(p1.id)
        pep3 = createMockPeptide(p1.id)
        pep4 = createMockPeptide(p1.id)
        mod = createMockPTM()
        
        createMockPeptideModification(m1, pep1, mod)
        createMockPeptideModification(m1, pep2, mod)
        createMockPeptideModification(m2, pep2, mod)
        createMockPeptideModification(m2, pep3, mod)
        createMockPeptideModification(m3, pep4, mod)
        
        measurements = [m1,m2,m3]
        
        pred1 = createMockScansite(p1.id, pep1.site_pos)
        pred2 = createMockScansite(p1.id, pep1.site_pos)
        pred3 = createMockScansite(p1.id, pep2.site_pos)
        pred4 = createMockScansite(p1.id, pep2.site_pos)
        pred5 = createMockScansite(p1.id, pep2.site_pos)
        pred6 = createMockScansite(p1.id, pep3.site_pos)
        
        pep1.predictions.extend([pred1, pred2])
        pep2.predictions = [pred3, pred4, pred5]
        pep3.predictions.extend([pred6])
        
        predictions = [pred1, pred2, pred3, pred4, pred5, pred6]
        
        pred1.source = 'scansite'
        pred2.source = 'scansite_bind'
        
        pred3.value = pred1.value
        pred6.value = pred1.value
        
        pred3.source = pred1.source
        pred6.source = pred1.source
        
        pred4.value = pred2.value
        pred4.source = pred2.source
        
        pred5.source = 'somewhere'
        
        return exp, measurements, predictions
    
    
    def test_format_predictions(self):
        _, measurements, predictions = self.createMockMeasurements(1)
        
        formatted_predictions = {}
        
        fstring = strings.prediction_type_map[predictions[0].source]         
        data = [(predictions[0].value, 2), ("None", 1)]
        jsondump = base64.b64encode(json.dumps(data))
        formatted_predictions[fstring] = {'json':jsondump, 'table':data}
        
        fstring = strings.prediction_type_map[predictions[1].source]
        data = [(predictions[1].value, 2), ("None", 1)]
        jsondump = base64.b64encode(json.dumps(data))
        formatted_predictions[fstring] = {'json':jsondump, 'table':data}
        
        data = [("None", 2), (predictions[4].value, 1)]
        jsondump = base64.b64encode(json.dumps(data))
        formatted_predictions[predictions[4].source] = {'json':jsondump, 'table':data}

        result = format_predictions(measurements)
        
        self.assertEqual(formatted_predictions, result)
    
    @patch('ptmscout.views.experiment.prediction_view.format_predictions')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByExperiment')
    @patch('ptmscout.database.experiment.getExperimentById')
    def test_prediction_view_should_format_data(self, patch_getExperiment, patch_getMeasurements, patch_format):
        request = DummyRequest()
        expid=1
        request.matchdict['id'] = expid
        request.user = createMockUser()
        
        exp, measurements, _ = self.createMockMeasurements(expid)
        
        patch_format.return_value = ["some formatted predictions"]
        patch_getExperiment.return_value = exp
        patch_getMeasurements.return_value = measurements
        
        request.user.experimentOwner.return_value = True
        
        result = prediction_view(request)
        
        request.user.experimentOwner.assert_called_once_with(exp)
        patch_getExperiment.assert_called_once_with(expid, request.user)
        patch_getMeasurements.assert_called_once_with(expid, request.user)
        patch_format.assert_called_once_with(measurements)
        
        
        self.assertEqual(strings.experiment_prediction_page_title % exp.name, result['pageTitle'])
        self.assertEqual(True, result['user_owner'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(patch_format.return_value, result['predictions'])

class TestPredictionViewIntegration(IntegrationTestCase):
    def test_prediction_view_integration(self):
        self.ptmscoutapp.get('/experiments/26/predictions')
