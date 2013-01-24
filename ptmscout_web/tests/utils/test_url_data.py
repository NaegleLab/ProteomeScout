import unittest
from pyramid.testing import DummyRequest
from ptmscout.utils.url_data import URLBuilder
import base64
import pickle
import urlparse

class TestURLBuilder(unittest.TestCase):

    def test_encode_url(self):
        request = DummyRequest()

        request.GET['experiment_id'] = '14'
        request.GET['site_pos'] = '234'
        request.GET['object'] = base64.urlsafe_b64encode( pickle.dumps( {"some":"kind","of":"thing"} ) )
        request.GET['accession'] = "ACK1_HUMAN"

        url_filter = URLBuilder()

        url_filter.create_field('experiment_id', URLBuilder.INTEGER)
        url_filter.create_field('object', URLBuilder.OBJECT)
        url_filter.create_field('accession', URLBuilder.STRING)
        url_filter.create_field('protein', URLBuilder.STRING)

        url_filter.decode_request(request)

        encoded_url = url_filter.encode_url( "/proteins/35546/modifications" )

        base, params = encoded_url.split("?")

        self.assertEqual('http://example.com/proteins/35546/modifications', base)
        decoded_params = urlparse.parse_qs(params)

        self.assertEqual({'experiment_id': [ '14' ], 'object':[ request.GET['object'] ], 'accession': ["ACK1_HUMAN"]},  decoded_params)

    def test_encode_url_empty(self):

        request = DummyRequest()
        url_filter = URLBuilder()
        url_filter.create_field('experiment_id', URLBuilder.INTEGER)
        
        url_filter.decode_request(request)
        encoded_url = url_filter.encode_url( "/proteins/35546/modifications" )

        self.assertEqual('http://example.com/proteins/35546/modifications', encoded_url)
