from mock import Mock
from ptmscout.utils import crypto
from ptmscout.database.user import User
from ptmscout.database.experiment import Experiment, ExperimentData
from ptmscout.database.permissions import Permission
from ptmscout.database.protein import Protein, GeneOntology, Domain
import random
from ptmscout.database.modifications import MeasuredPeptide, Phosphopep,\
    ScansitePrediction
from ptmscout.database.gene_expression import ExpressionProbeset,\
    ExpressionSample, ExpressionCollection, ExpressionTissue
from ptmscout.database.taxonomies import Species
from ptmscout.database.upload import Session, SessionColumn
import datetime

TEST_USER_ID = 2

def createMockUser(username=None, email=None, password=None, active=1):
    global TEST_USER_ID
    mock = Mock(spec=User)
    
    mock.id = TEST_USER_ID
    
    if username==None:
        username = "user" + str(TEST_USER_ID)
    mock.username = username
    mock.name = "A User"
    mock.email = email
    
    if email == None:
        email = "user%d@institute.edu" % TEST_USER_ID
    mock.email = email
    
    mock.institution = "institution"
    
    TEST_USER_ID += 1
    if password==None:
        password=crypto.randomString(10)
    
    mock.salt, mock.salted_password = crypto.saltedPassword(password)  
    mock.activation_token = crypto.generateActivationToken()
    mock.active = active
    mock.permissions = []
    return mock

def createMockSession(user, sid=random.randint(0,100000), data_file='some_file', load_type='new', parent_experiment=None, stage = 'config', experiment_id=None, change_description='blah'):
    mock = Mock(spec=Session)
    mock.id=sid
    mock.user_id=user.id
    mock.data_file = data_file
    mock.load_type = load_type
    mock.parent_experiment = parent_experiment
    mock.change_description = change_description
    mock.stage = stage
    mock.experiment_id = experiment_id
    mock.columns = []
    mock.date = datetime.datetime.now()
    return mock
    
def createMockSessionColumn(col_num, tp, session_id, label='', scid=random.randint(0,100000)):
    mock = Mock(spec=SessionColumn)
    mock.id = scid
    mock.session_id = session_id
    mock.type = tp
    mock.label = label
    mock.column_number = col_num
    return mock

def createMockExperiment(eid=random.randint(0,100000), public=0, parent_id=None, status='loaded'):
    mock = Mock(spec=Experiment)
    mock.id = eid
    mock.public = public
    mock.parent_id = parent_id
    mock.URL = "url"
    mock.getUrl.return_value = mock.URL
    mock.name = "Experiment Name" + str(eid)
    mock.status = status
    mock.ready.return_value = mock.status == 'loaded'
    mock.date = datetime.datetime.now()
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
    
    pid = random.randint(0,100000)
    mock.id = pid
    mock.name = "prot_" + str(pid)
    mock.acc_gene = "PR" + str(pid)
    mock.date="12-1986"
    mock.species_id=46
    mock.species = mock(spec=Species)
    mock.species.id = 46
    mock.species.name = "homo sapiens"
    mock.sequence="ABCDEFGHIJKLMNOP" 
    mock.GO_terms = []
    return mock

def createMockGO(go_type):
    mock = Mock(spec=GeneOntology)
    
    gid = random.randint(0,100000)
    mock.id = gid
    
    if(go_type not in ['F','P','C']):
        raise ValueError("Bad GO type: " + str(go_type))
    
    mock.aspect = go_type
    
    mock.GO = "GO:" + str(gid)
    mock.term = "Some term " + str(gid)
    mock.version = "1.2"
    
    return mock
    
def createMockMeasurement(pid, expid):
    mock = Mock(spec=MeasuredPeptide)
    
    mock.id = random.randint(0,100000)
    mock.protein_id = pid
    mock.experiment_id = expid
    mock.phosphopep = "ABCDEF"
    mock.phosphopeps = []
    mock.data = []
    
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
    
    mock.predictions = []
    
    return mock

def createMockData(num, run, mod_id):
    rval = []
    
    labels = [0,1,5,10,20,50,100,200,500]
    
    for i in xrange(0, num):
        mock = Mock(spec=ExperimentData)
        
        mock.id = random.randint(0, 100000)
        
        mock.type = 'time'
        mock.run = run
        mock.label = labels[i] 
        
        mock.priority = i+1
        mock.value = random.random();
        
        mock.MS_id = mod_id
        
        rval.append(mock)
    
    random.shuffle(rval)
    return rval

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
    
def createMockScansite(pep_id):
    mock = Mock(spec=ScansitePrediction)
    
    mock.id = random.randint(0, 100000)
    mock.source = "somesource" + str(mock.id)
    mock.value = "Some_Value" + str(mock.id)
    mock.score = random.random()
    mock.phophopep_id = pep_id
    
    return mock

def createMockDomain(pid):
    mock = Mock(spec=Domain)
    
    mock.id = random.randint(0, 100000)
    
    mock.label = "some_label"+str(mock.id)
    mock.start = random.randint(0, 10000)
    mock.stop = random.randint(mock.start+1, 10000)
    
    mock.p_value = random.random()
    
    mock.source = 'pfam'
    mock.params = 'pval=1e-05'
    
    mock.protein_id = pid
    mock.version = 23
    
    return mock
     
    
    
    