from sqlalchemy import Column, Integer, Text

class PTMUser(object):
    __tablename__ = 'users'
    id = Column(Integer, primary_key=True)
    username = Column(Text, unique=True)
    salted_password = Column(Text)
    name = Column(Text)
    email = Column(Text, unique=True)
    institution = Column(Text)
    
    def init(self, username, name, email, institution):
        pass
    
    def createUser(self, password):
        pass
