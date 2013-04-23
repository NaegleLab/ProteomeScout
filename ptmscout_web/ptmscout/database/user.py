from . import Base, DBSession
from sqlalchemy import Column, Integer, VARCHAR, TIMESTAMP
import ptmscout.utils.crypto as crypto 
from pyramid import security
from sqlalchemy.orm import relationship
from ptmscout.database import permissions
from sqlalchemy.types import Enum, DateTime
import datetime

class User(Base):
    __tablename__ = 'users'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    # credentials
    username = Column(VARCHAR(30), unique=True)
    salted_password = Column(VARCHAR(64))
    salt = Column(VARCHAR(10))
    
    # user data
    name = Column(VARCHAR(50))
    email = Column(VARCHAR(50), unique=True)
    institution = Column(VARCHAR(100))
    date_created = Column(TIMESTAMP)
    
    # activation data
    active = Column(Integer(1), default=0)
    activation_token = Column(VARCHAR(50))
    
    access_level = Column(Enum(['reviewer','researcher']), default='researcher')
    expiration = Column(DateTime)

    permissions = relationship("Permission", backref="user")
    
    def __init__(self, username="", name="", email="", institution="", access_level='researcher'):
        self.username = username
        self.name = name
        self.email = email
        self.institution = institution
        self.access_level = access_level
        self.permissions = []
        
    def isExpired(self):
        if self.expiration != None:
            now = datetime.datetime.now()
            if self.expiration < now:
                return True
        return False    
        
    def setExpiration(self, delta):
        n = datetime.datetime.now()
        d = datetime.timedelta(delta)
        self.expiration = n + d
        
    def saveUser(self):
        DBSession.add(self)
        DBSession.flush()
        
    def setActive(self):
        self.active = 1
        
    def createUser(self, password):
        self.salt, self.salted_password = crypto.saltedPassword(password)
        self.activation_token = crypto.generateActivationToken()
        
    def myExperiments(self):
        return [ permission.experiment \
                    for permission in self.permissions \
                        if permission.access_level == 'owner' and permission.experiment.isExperiment() ]
        
    def myDatasets(self):
        return [ permission.experiment \
                    for permission in self.permissions \
                        if permission.access_level == 'owner' and not permission.experiment.isExperiment() ]

        
    def experimentOwner(self, exp):
        result = [ permission.experiment \
                        for permission in self.permissions \
                            if permission.access_level == 'owner' and permission.experiment.id == exp.id]
        return len(result) > 0
    
    def allExperiments(self):
        return [ permission.experiment for permission in self.permissions ]

    
    def processInvitations(self):
        invites = permissions.getInvitationsForUser(self.email)
        
        for invite in invites:
            np = permissions.Permission(invite.experiment)
            np.user_id = self.id
            self.permissions.append(np)
        
    
    
        
class NoSuchUser(Exception):
    def __init__(self, uid=None, username=None, email=None):
        self.uid = uid
        self.username = username
        self.email = email
    def __str__(self):
        value = ""
        if self.uid != None:
            value = str(self.uid)
        if self.username != None:
            value = str(self.username)
        if self.email != None:
            value = str(self.email)
        
        return "Could not find user %s" % value

def getUserByRequest(request):
    username = security.authenticated_userid(request)
    if username is not None:
        try:
            return getUserByUsername(username)
        except NoSuchUser:
            return None
    return None

def getUsernameByRequest(request):
    return security.authenticated_userid(request)

def getUserById(uid):
    value = DBSession.query(User).filter_by(id=uid).first()
    if value == None:
        raise NoSuchUser(uid=uid)
    return value
    
def getUserByUsername(username):
    value = DBSession.query(User).filter_by(username=username).first()
    if value == None:
        raise NoSuchUser(username=username)
    return value

def getUserByEmail(email):
    value = DBSession.query(User).filter_by(email=email).first()
    if value == None:
        raise NoSuchUser(email=email)
    return value
    
    
