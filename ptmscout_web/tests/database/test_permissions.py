from tests.DBTestCase import DBTestCase
import ptmscout.database.user as user
import ptmscout.database.experiment as experiment
from ptmscout.database import permissions

class PermissionsTestCase(DBTestCase):
    def test_permission_should_get_correct_permissions(self):
        u = user.User("user","name","email","inst")
        u.createUser("password")
        
        e = experiment.getExperimentById(1, None)
        
        p = permissions.Permission(e, 'owner')
        
        u.permissions.append(p)
        u.saveUser()
        
        self.assertEqual(u.id, experiment.getExperimentById(1, u).permissions[0].user.id)
        
    def test_getInvitationsForUser_should_return_correct_invites(self):
        inviter = user.User("user","name","email","inst")
        inviter.createUser("password")
        inviter.saveUser() 
        
        i1 = permissions.Invitation("invited_user@institute.edu", 1, inviter.id)
        i2 = permissions.Invitation("invited_user@institute.edu", 26, inviter.id)
        i3 = permissions.Invitation("invited_user@institute.edu", 28, inviter.id)
        
        [ invite.saveInvitation() for invite in [i1,i2,i3] ]
        
        invitation_list = permissions.getInvitationsForUser("invited_user@institute.edu")
        invited_exp_ids = [ invite.experiment_id for invite in invitation_list ]
        self.assertEqual([1,26,28], invited_exp_ids)
