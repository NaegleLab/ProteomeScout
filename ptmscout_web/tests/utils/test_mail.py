import unittest
from ptmscout.utils import mail

class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def test_email_is_valid_valid_email_not_valid_domain(self):
        ok, domain_ok  = mail.email_is_valid('user@example.com')
        
        self.assertEqual(True, ok)
        self.assertEqual(False, domain_ok)
    
    def test_email_is_valid_not_valid_email(self):
        ok, domain_ok  = mail.email_is_valid('user@example')
        
        self.assertEqual(False, ok)
        self.assertEqual(False, domain_ok)
    
    def test_email_is_valid_valid(self):
        ok, domain_ok  = mail.email_is_valid('user@example.edu')
        
        self.assertEqual(True, ok)
        self.assertEqual(True, domain_ok)
        
        
        
