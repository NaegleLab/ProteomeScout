from models import Base
from sqlalchemy import Column, Integer, VARCHAR, TIMESTAMP
import hashlib
import os
from sqlalchemy.util import buffer

class PTMUser(Base):
    __tablename__ = 'users'
    id = Column(Integer(10), primary_key=True)
    
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
    active = Column(Integer(1), default = 0)
    activation_token = Column(VARCHAR(50))
    
    def __init__(self, username, name, email, institution):
        self.username = username
        self.name = name
        self.email = email
        self.institution = institution
        
    def createUser(self, password):
        self.salt = self.__byteStringToHex(os.urandom(5))
        self.salted_password = hashlib.sha256(self.salt + password).hexdigest()
        self.activation_token = self.__byteStringToHex(os.urandom(25))
        
    def checkPassword(self, password):
        pass
    
    def __byteStringToHex(self, byteString):
        return ''.join([hex(ord(c))[2:].zfill(2) for c in byteString])