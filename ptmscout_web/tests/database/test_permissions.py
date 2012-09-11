from tests.DBTestCase import DBTestCase
import ptmscout.database.user as user
import ptmscout.database.experiment as experiment
from ptmscout.database import permissions

class PermissionsTestCase(DBTestCase):
    def test_permission_should_get_correct_permissions(self):
        u = user.User("user","name","email","inst")
        u.createUser("password")
        
        e = experiment.getExperimentById(1)
        
        p = permissions.Permission(e, 'owner')
        
        u.permissions.append(p)
        u.saveUser()
        
        self.assertEqual(u.id, experiment.getExperimentById(1).permissions[0].user.id)