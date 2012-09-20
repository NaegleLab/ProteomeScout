from mock import Mock
from ptmscout.utils import crypto
from ptmscout.database.user import User
from ptmscout.database.experiment import Experiment
from ptmscout.database.permissions import Permission
from ptmscout.database.protein import Protein, Species
import random
from ptmscout.database.modifications import Modification, Phosphopep
from ptmscout.database.gene_expression import ExpressionProbeset,\
    ExpressionSample, ExpressionCollection, ExpressionTissue

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
    mock.species_id=46
    mock.species = mock(spec=Species)
    mock.species.id = 46
    mock.species.name = "homo sapiens"
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
    
    mock.getName.return_value = mock.site_type + str(mock.site_pos)
    mock.getPeptide.return_value = mock.pep_aligned
    
    return mock

def createMockProbe():
    mock = Mock(spec=ExpressionProbeset)
    
    mock.id = random.randint(0, 100000)
    mock.probeset_id = str(mock.id) + "_s_at"
    st = ['gnf1h', 'gnf1m', 'HG-U133A']
    mock.genechip = st[random.randint(0, 2)]
    
    mock.species = mock(spec=Species)
    mock.species.id = 46
    mock.species.name = "homo sapiens"
    
    mock.name = "Some probeset"
    mock.samples = []
    mock.accessions = []
    return mock
    
def createMockExpSample(probe_id, collection_id, collection, tissue_id, tissue):
    mock = Mock(spec=ExpressionSample)
    
    mock.id = random.randint(0, 100000)
    mock.probeset_id = probe_id

    mock.collection_id = collection_id
    mock.collection = Mock(spec=ExpressionCollection)
    mock.collection.id = collection_id
    mock.collection.name = collection
    
    mock.tissue_id = tissue_id
    mock.tissue = Mock(spec=ExpressionTissue)
    mock.tissue.id = tissue_id
    mock.tissue.name = tissue
    
    mock.value = random.random()
    return mock
    