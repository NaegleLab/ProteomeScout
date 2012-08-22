from pyramid import testing
from pyramid.httpexceptions import HTTPFound
from pyramid.testing import DummyRequest
import unittest
import urllib

class UserManagementTests(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()

    def test_login_should_display_login_page(self):
        from user_management import user_login
        request = DummyRequest()
    
        info = user_login(request)
        
        self.assertEqual('Login', info['pageTitle'])
        self.assertEqual(None, info['reason'])
        
    def test_login_should_display_error_with_reason_and_populate_username(self):
        from user_management import user_login
        request = DummyRequest()
        
        request.GET['username'] = "a_username"
        request.GET['reason'] = "you didn't do it right"
        
        info = user_login(request)
        
        self.assertEqual('Login', info['pageTitle'])
        self.assertEqual('a_username', info['username'])
        self.assertEqual("you didn't do it right", info['reason'])
        
    def test_process_login_should_redirect_to_login_with_error_on_fields_missing(self):
        from user_management import user_login_success
        request = DummyRequest()
        
        request.POST['username'] = "a_username"
        request.POST['password'] = ""
        
        try:
            user_login_success(request)
            self.fail("Expected exception HTTPFound, no exception raised")
        except HTTPFound, f:
            self.assertEqual("http://example.com/login?"+urllib.urlencode({'username':"a_username",'reason':"All fields are required"}), f.location)
            
    def test_process_login_should_fail_if_credentials_incorrect(self):
        from user_management import user_login_success
        request = DummyRequest()
        
        request.GET['username'] = "good_username"
        request.GET['password'] = "bad_password"

        try:
            user_login_success(request)
            self.fail("Expected exception HTTPFound, no exception raised")
        except HTTPFound, f:
            self.assertEqual("http://example.com/login?"+urllib.urlencode({'username':"a_username",'reason':"Credentials incorrect"}), f.location)
    
    def test_process_login_should_fail_if_no_such_user(self):
        from user_management import user_login_success
        request = DummyRequest()
        
        request.GET['username'] = "notauser"
        request.GET['password'] = "password"
        
        try:
            user_login_success(request)
            self.fail("Expected exception HTTPFound, no exception raised")
        except HTTPFound, f:
            self.assertEqual("http://example.com/login?"+urllib.urlencode({'username':"a_username",'reason':"Credentials incorrect"}), f.location)
        