from tests.DBTestCase import DBTestCase
import database.user
from database.user import PTMUser
from sqlalchemy.util import buffer


class UserTestCase(DBTestCase):
    def test_createUser_should_insert_user(self):
        user = PTMUser("newguy", "new guy's name", "newguy@someschool.edu", "some school")
        user.createUser("password")
        user.saveUser()
        
        newuser = database.user.getUserByUsername("newguy")
        self.assertEqual(user.username, newuser.username)
        
        self.assertEqual(user.salted_password, newuser.salted_password)
        self.assertEqual(user.salt, newuser.salt)
        
        self.assertEqual(user.name, newuser.name)
        self.assertEqual(user.email, newuser.email)
        self.assertEqual(user.institution, newuser.institution)
        self.assertEqual(user.activation_token, newuser.activation_token)
        
        self.assertEqual(0, user.active)
        self.assertNotEqual(None, newuser.date_created)
        
    def test_createUser_should_set_correct_id(self):
        user = PTMUser("newguy", "new guy's name", "newguy@someschool.edu", "some school")
        user.createUser("password")
        user.saveUser()
        
        newuser = database.user.getUserById(user.id)
        
        self.assertEqual(user.username, newuser.username) 