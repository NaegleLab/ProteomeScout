from tests.DBTestCase import DBTestCase
import ptmscout.database.user as dbuser
from ptmscout.database.user import User, NoSuchUser, getUserById
from ptmscout.database.experiment import Experiment


class UserTestCase(DBTestCase):
    def test_user_should_associate_with_correct_permissions(self):
        experiments = [Experiment(),Experiment(),Experiment()] 
        for i, exp in enumerate(experiments):
            exp.id = 1000 + i
            exp.name = "experiment " + str(exp.id)
            exp.author = ""
            exp.date = ""
            exp.published = ""
            exp.ambiguity = 0
            exp.export = 0
            exp.experiment_id = 0
            exp.dataset = ""
            exp.submitter = ""
            exp.primaryModification = "Y"
            exp.saveExperiment()
            
        users = [User(), User()]
        
        for i, user in enumerate(users):
            user.id = 100+i
            user.username = "user" + str(i)
            user.createUser("password")
            user.active=1
        
        users[0].experiments.append(experiments[0])
        users[0].experiments.append(experiments[1])
        
        users[1].experiments.append(experiments[1])
        users[1].experiments.append(experiments[2])
        
        [user.saveUser() for user in users]
        
        self.assertEqual([experiments[0].id, experiments[1].id], [exp.id for exp in getUserById(100).experiments])
        self.assertEqual([experiments[1].id, experiments[2].id], [exp.id for exp in getUserById(101).experiments])
        
        
        
    def test_createUser_should_insert_user(self):
        user = User("newguy", "new guy's name", "newguy@someschool.edu", "some school")
        user.createUser("password")
        user.saveUser()

        try:        
            newuser = dbuser.getUserById(user.id)
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
        user = User("newguy", "new guy's name", "newguy@someschool.edu", "some school")
        user.createUser("password")
        user.saveUser()
        
        newuser = dbuser.getUserById(user.id)
        
        self.assertEqual(user.username, newuser.username) 
        
    def test_getUserByUsername_should_throw_exception_when_no_such_user(self):
        try:
            dbuser.getUserByUsername("nosuchuser")
        except NoSuchUser, n:
            self.assertEqual("nosuchuser", n.username)
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception NoSuchUser was not thrown")
        