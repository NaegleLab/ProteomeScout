from tests.DBTestCase import DBTestCase
import ptmscout.database.user as dbuser
from ptmscout.database.user import PTMUser, NoSuchUser


class UserTestCase(DBTestCase):
    def test_createUser_should_insert_user(self):
        user = PTMUser("newguy", "new guy's name", "newguy@someschool.edu", "some school")
        user.createUser("password")
        user.saveUser()

        try:        
            newuser = dbuser.getUserByUsername("newguy")
        except NoSuchUser:
            self.fail("User insert failed")
            
        self.assertEqual(user.username, newuser.username)
        
        self.assertEqual(user.salted_password, newuser.salted_password)
        self.assertEqual(user.salt, newuser.salt)
        
        self.assertEqual(user.name, newuser.name)
        self.assertEqual(user.email, newuser.email)
        self.assertEqual(user.institution, newuser.institution)
        self.assertEqual(user.activation_token, newuser.activation_token)
        
        self.assertEqual(0, user.active)
        self.assertEqual(0, newuser.active)
        
        self.assertNotEqual(None, user.date_created)
        self.assertEqual(user.date_created, newuser.date_created)
        
    def test_createUser_should_set_correct_id(self):
        user = PTMUser("newguy", "new guy's name", "newguy@someschool.edu", "some school")
        user.createUser("password")
        user.saveUser()
        
        newuser = dbuser.getUserById(user.id)
        
        self.assertEqual(user.username, newuser.username) 