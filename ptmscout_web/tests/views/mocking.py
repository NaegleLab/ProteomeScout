from mock import Mock
from ptmscout.utils import crypto
from ptmscout.database.user import User
from ptmscout.database.experiment import Experiment
from ptmscout.database.permissions import Permission

TEST_USER_ID = 17

def createUserForTest(username, email, password, active):
    mock = Mock(spec=User)
    mock.username = username
    mock.name = "A User"
    mock.email = email
    mock.institution = "institution"
    mock.salt, mock.salted_password = crypto.saltedPassword(password)  
    mock.activation_token = crypto.generateActivationToken()
    mock.id = TEST_USER_ID
    mock.active = active
    return mock

def createMockExperiment(eid, public, parent_id=0):
    mock = Mock(spec=Experiment)
    mock.id = eid
    mock.public = public
    mock.parent_id = parent_id
    return mock

def createMockPermission(user, experiment, access_level='view'):
    mock = Mock(spec=Permission)
    
    mock.user = user
    mock.experiment = experiment
    mock.user_id = user.id
    mock.experiment_id = experiment.id
    mock.access_level = access_level
    return mock