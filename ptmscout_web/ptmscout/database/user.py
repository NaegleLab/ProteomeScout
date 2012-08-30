from . import Base, DBSession
from sqlalchemy import Column, Integer, VARCHAR, TIMESTAMP
import ptmscout.utils.crypto as crypto 
from pyramid import security

class PTMUser(Base):
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
    
    def __init__(self, username, name, email, institution):
        self.username = username
        self.name = name
        self.email = email
        self.institution = institution
        
    def saveUser(self):
        DBSession.add(self)
        DBSession.flush()
        
    def createUser(self, password):
        self.salt, self.salted_password = crypto.saltedPassword(password)  
        self.activation_token = crypto.generateActivationToken()
        
    def checkPassword(self, password):
        pass
    
    
class NoSuchUser(Exception):
    def __init__(self, uid=None, username=None):
        self.uid = uid
        self.username = username
    def __str__(self):
        value = ""
        if self.uid != None:
            value = str(self.uid)
        if self.username != None:
            value = str(self.username)
        
        return "Could not find user %s" % value

def getUserByRequest(request):
    username = security.authenticated_userid(request)
    if username is not None:
        try:
            return getUserByUsername(username)
        except NoSuchUser:
            return None

def getUserById(uid):
    value = DBSession.query(PTMUser).filter_by(id=uid).first()
    if value == None:
        raise NoSuchUser(uid=uid)
    return value
    
def getUserByUsername(username):
    value = DBSession.query(PTMUser).filter_by(username=username).first()
    if value == None:
        raise NoSuchUser(username=username)
    return value
    
    
