from . import Base, DBSession
from sqlalchemy import Column, Integer, VARCHAR, TIMESTAMP
import ptmscout.utils.crypto as crypto 
from pyramid import security
from sqlalchemy.orm import relationship


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

    permissions = relationship("Permission", backref="user")
    
    def __init__(self, username="", name="", email="", institution=""):
        self.username = username
        self.name = name
        self.email = email
        self.institution = institution
        
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
                        if permission.access_level == 'owner' ]
        
    def allExperiments(self):
        return [ permission.experiment for permission in self.permissions ]
        
    def makeOwner(self, exp_id):
        DBSession.execute("UPDATE permissions SET access_level='owner' WHERE user_id=%d and experiment_id=%d" % (self.id, exp_id))        
        
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
    
    
