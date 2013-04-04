from tests.DBTestCase import DBTestCase
import ptmscout.database.user as dbuser
from ptmscout.database.user import User, NoSuchUser, getUserById
from ptmscout.database.experiment import Experiment
from ptmscout.database.permissions import Permission
from ptmscout.database import permissions


class UserTestCase(DBTestCase):
    def test_user_should_associate_with_correct_permissions(self):
        experiments = [Experiment(),Experiment(),Experiment()] 
        for i, exp in enumerate(experiments):
            exp.id = 1000 + i
            exp.name = "experiment " + str(exp.id)
            exp.author = ""
            exp.published = 0
            exp.ambiguity = 0
            exp.type = 'compendia'
            exp.experiment_id = None
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
        
        
        users[0].permissions.append(Permission(experiments[0]))
        users[0].permissions.append(Permission(experiments[1]))
        
        users[1].permissions.append(Permission(experiments[1]))
        users[1].permissions.append(Permission(experiments[2]))
        
        [user.saveUser() for user in users]
        
        self.assertEqual([experiments[0].id, experiments[1].id], [p.experiment.id for p in getUserById(100).permissions])
        self.assertEqual([experiments[1].id, experiments[2].id], [p.experiment.id for p in getUserById(101).permissions])
        
        
        
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
            
    def test_getUserByEmail_should_return_user_when_exists(self):
        user = User("newguy", "new guy's name", "newguy@someschool.edu", "some school")
        user.createUser("password")
        user.saveUser()
        
        newuser = dbuser.getUserByEmail("newguy@someschool.edu")
        
        self.assertEqual(user.username, newuser.username) 

    def test_getUserByEmail_should_throw_exception_when_no_such_user(self):
        try:
            dbuser.getUserByEmail("nosuchuser")
        except NoSuchUser, n:
            self.assertEqual("nosuchuser", n.email)
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception NoSuchUser was not thrown")
            
    def test_processInvitations_should_query_for_invites_and_add_experiments_to_permissions(self):
        inviter = User("inviter", "inviter's name", "inviter@someschool.edu", "some school")
        inviter.createUser("password")
        inviter.saveUser()
        
        new_user = User("newguy", "new guy's name", "newguy@someschool.edu", "some school")
        new_user.createUser("password")
        new_user.saveUser()
        
        invites = [permissions.Invitation(new_user.email, 26, inviter.id), 
                   permissions.Invitation(new_user.email, 28, inviter.id)]
        
        [ i.saveInvitation() for i in invites ]
        self.session.flush()
        
        new_user.processInvitations()
        new_user.saveUser()
        
        user = getUserById(new_user.id)
        
        allowed_exp_ids = [ p.experiment.id for p in user.permissions ]
        
        self.assertEqual([26,28], allowed_exp_ids)
