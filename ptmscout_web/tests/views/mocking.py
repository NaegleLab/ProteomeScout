from mock import Mock
from ptmscout.utils import crypto
from ptmscout.database.user import User
from ptmscout.database.experiment import Experiment
from ptmscout.database.permissions import Permission
from ptmscout.database.protein import Protein
import random
from ptmscout.database.modifications import Modification, Phosphopep

TEST_USER_ID = 2

def createUserForTest(username, email, password, active):
    global TEST_USER_ID
    mock = Mock(spec=User)
    mock.username = username
    mock.name = "A User"
    mock.email = email
    mock.institution = "institution"
    mock.salt, mock.salted_password = crypto.saltedPassword(password)  
    mock.activation_token = crypto.generateActivationToken()
    mock.id = TEST_USER_ID
    TEST_USER_ID += 1
    mock.active = active
    mock.permissions = []
    return mock

def createMockExperiment(eid, public, parent_id=0):
    mock = Mock(spec=Experiment)
    mock.id = eid
    mock.public = public
    mock.parent_id = parent_id
    mock.name = "Experiment Name"
    return mock

def createMockPermission(user, experiment, access_level='view'):
    mock = Mock(spec=Permission)
    
    mock.user = user
    mock.experiment = experiment
    mock.user_id = user.id
    mock.experiment_id = experiment.id
    mock.access_level = access_level
    return mock

def createMockProtein():
    mock = Mock(spec=Protein)
    
    id = random.randint(0,100000)
    mock.id = id
    mock.name = "prot_" + str(id)
    mock.acc_gene = "PR" + str(id)
    mock.date="12-1986"
    mock.species="homo sapiens"
    mock.sequence="ABCDEFGHIJKLMNOP" 
    return mock
    
def createMockModification(pid, expid):
    mock = Mock(spec=Modification)
    
    mock.id = random.randint(0,100000)
    mock.protein_id = pid
    mock.experiment_id = expid
    mock.phosphopep = "ABCDEF"
    mock.phosphopeps = []
    
    return mock

def createMockPhosphopep(pid):
    mock = Mock(spec=Phosphopep)
    
    site_types = ['Y','T','S','P']
    
    mock.id = random.randint(0, 100000)
    mock.protein_id = pid
    mock.pep_aligned = "ABCDEFHG"
    mock.pep_tryps = "JFDAABCDEFHGADFG"
    mock.site_pos = random.randint(0, 3000)
    mock.site_type = site_types[random.randint(0,3)]
    mock.pfam_site = "pfam1"
    
    return mock