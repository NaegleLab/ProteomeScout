import unittest
from pyramid import testing
from pyramid.testing import DummyRequest
from tests.views.mocking import createMockProtein, createMockMeasurement,\
    createMockPhosphopep, createMockExperiment, createMockScansite
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
        pred1 = createMockScansite(1)
        pred2 = createMockScansite(1)
        pred3 = createMockScansite(1)
        
        pred1.score = 3
        pred2.score = 0.2
        pred3.score = 0.1
        
        pred3.value = "~~~"
        
        result = filter_predictions([pred1,pred2,pred3])
        
        self.assertEqual([pred2], result)
        

    def createMockMeasurements(self, expid):
        exp = createMockExperiment(expid, 0, 0)
        
        p1 = createMockProtein()
        m1 = createMockMeasurement(p1.id, expid)
        m2 = createMockMeasurement(p1.id, expid)
        
        pep1 = createMockPhosphopep(p1.id)
        pep2 = createMockPhosphopep(p1.id)
        pep3 = createMockPhosphopep(p1.id)
        
        m1.phosphopeps.extend([pep1, pep2])
        m2.phosphopeps.extend([pep2, pep3])
        
        measurements = [m1,m2]
        
        pred1 = createMockScansite(pep1.id)
        pred2 = createMockScansite(pep1.id)
        pred3 = createMockScansite(pep2.id)
        pred4 = createMockScansite(pep2.id)
        pred5 = createMockScansite(pep2.id)
        pred6 = createMockScansite(pep3.id)
        
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
        data = [(predictions[0].value, 2)]
        jsondump = base64.b64encode(json.dumps(data))
        formatted_predictions[fstring] = {'json':jsondump, 'table':data}
        
        fstring = strings.prediction_type_map[predictions[1].source]
        data = [(predictions[1].value, 2)]
        jsondump = base64.b64encode(json.dumps(data))
        formatted_predictions[fstring] = {'json':jsondump, 'table':data}
        
        data = [(predictions[4].value, 1)]
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
        request.user = None
        
        exp, measurements, _ = self.createMockMeasurements(expid)
        
        patch_format.return_value = ["some formatted predictions"]
        patch_getExperiment.return_value = exp
        patch_getMeasurements.return_value = measurements
        
        result = prediction_view(request)
        
        patch_getExperiment.assert_called_once_with(expid, request.user)
        patch_getMeasurements.assert_called_once_with(expid, request.user)
        patch_format.assert_called_once_with(measurements)
        
        
        self.assertEqual(strings.experiment_prediction_page_title % exp.name, result['pageTitle'])
        self.assertEqual(exp, result['experiment'])
        self.assertEqual(patch_format.return_value, result['predictions'])

class TestPredictionViewIntegration(IntegrationTestCase):
    def test_prediction_view_integration(self):
        self.ptmscoutapp.get('/experiments/26/predictions')
