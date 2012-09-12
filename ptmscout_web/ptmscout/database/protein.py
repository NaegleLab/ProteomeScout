from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column
from sqlalchemy.types import Integer, TEXT, VARCHAR


class Protein(Base):
    __tablename__='protein'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    sequence = Column(TEXT)
    species = Column(VARCHAR(100))
    acc_gene = Column(VARCHAR(30))
    name = Column(VARCHAR(100))
    date = Column(VARCHAR(7))
    
    def __init__(self):
        pass
    
    def saveProtein(self):
        DBSession.add(self)

class NoSuchProtein(Exception):
    def __init__(self, pid):
        self.pid = pid
    
    def __str__(self):
        return "No such protein: %d" % (self.pid)

def getProteinById(pid):
    value = DBSession.query(Protein).filter_by(id=pid).first()
    
    if value == None:
        raise NoSuchProtein(pid)
    
    return value